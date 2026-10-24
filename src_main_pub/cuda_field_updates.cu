#include <cuda_runtime.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

namespace {

enum FieldOperation {
    advance_ex = 1,
    advance_ey = 2,
    advance_ez = 3,
    advance_hx = 4,
    advance_hy = 5,
    advance_hz = 6
};

struct FieldView {
    void* data;
    int nx;
    int ny;
    int nz;
    int lo_x;
    int lo_y;
    int lo_z;
};

struct VectorView {
    void* data;
    int length;
    int lo;
};

struct Sweep {
    int nx;
    int ny;
    int nz;
};

struct CpmlJob {
    FieldView output;
    FieldView media;
    FieldView source;
    FieldView psi;
    VectorView coefficient_b;
    VectorView coefficient_c;
    VectorView coefficient_g;
    int region_xi;
    int region_xe;
    int region_yi;
    int region_ye;
    int region_zi;
    int region_ze;
    int axis;
    int source_direction;
    int update_sign;
};

class DeviceBuffer {
public:
    DeviceBuffer() = default;
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;

    ~DeviceBuffer() {
        if (data_ != nullptr) {
            cudaFree(data_);
        }
    }

    cudaError_t copy_from(const void* host_data, std::size_t bytes) {
        if (bytes == 0) {
            return cudaErrorInvalidValue;
        }
        auto status = cudaMalloc(&data_, bytes);
        if (status != cudaSuccess) {
            return status;
        }
        return cudaMemcpy(data_, host_data, bytes, cudaMemcpyHostToDevice);
    }

    cudaError_t copy_to(void* host_data, std::size_t bytes) const {
        return cudaMemcpy(host_data, data_, bytes, cudaMemcpyDeviceToHost);
    }

    template <typename T>
    T* as() const {
        return static_cast<T*>(data_);
    }

private:
    void* data_ = nullptr;
};

__host__ __device__ std::size_t index3(int i, int j, int k, int nx, int ny) {
    return static_cast<std::size_t>(i)
        + static_cast<std::size_t>(nx)
            * (static_cast<std::size_t>(j) + static_cast<std::size_t>(ny) * k);
}

template <typename Real>
__global__ void update_field_kernel(
        int operation,
        Real* output,
        const std::int16_t* media,
        const Real* field_a,
        const Real* field_b,
        const Real* metric_a,
        const Real* metric_b,
        const Real* coefficient_a,
        const Real* coefficient_b,
        int output_nx,
        int output_ny,
        int media_nx,
        int media_ny,
        int field_a_nx,
        int field_a_ny,
        int field_b_nx,
        int field_b_ny,
        int sweep_nx,
        int sweep_ny,
        int sweep_nz,
        int total) {
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total) {
        return;
    }

    const int i = 1 + index % sweep_nx;
    const int j = 1 + (index / sweep_nx) % sweep_ny;
    const int k = 1 + index / (sweep_nx * sweep_ny);
    const int medium = media[index3(i, j, k, media_nx, media_ny)];
    const int output_index = index3(i, j, k, output_nx, output_ny);
    const Real old_value = output[output_index];
    Real curl;

    switch (operation) {
    case advance_ex:
        curl = (field_a[index3(i, j, k, field_a_nx, field_a_ny)]
                - field_a[index3(i, j - 1, k, field_a_nx, field_a_ny)]) * metric_a[j]
            - (field_b[index3(i, j, k, field_b_nx, field_b_ny)]
                - field_b[index3(i, j, k - 1, field_b_nx, field_b_ny)]) * metric_b[k];
        break;
    case advance_ey:
        curl = (field_a[index3(i, j, k, field_a_nx, field_a_ny)]
                - field_a[index3(i, j, k - 1, field_a_nx, field_a_ny)]) * metric_a[k]
            - (field_b[index3(i, j, k, field_b_nx, field_b_ny)]
                - field_b[index3(i - 1, j, k, field_b_nx, field_b_ny)]) * metric_b[i];
        break;
    case advance_ez:
        curl = (field_a[index3(i, j, k, field_a_nx, field_a_ny)]
                - field_a[index3(i - 1, j, k, field_a_nx, field_a_ny)]) * metric_a[i]
            - (field_b[index3(i, j, k, field_b_nx, field_b_ny)]
                - field_b[index3(i, j - 1, k, field_b_nx, field_b_ny)]) * metric_b[j];
        break;
    case advance_hx:
        curl = (field_a[index3(i, j, k + 1, field_a_nx, field_a_ny)]
                - field_a[index3(i, j, k, field_a_nx, field_a_ny)]) * metric_a[k]
            - (field_b[index3(i, j + 1, k, field_b_nx, field_b_ny)]
                - field_b[index3(i, j, k, field_b_nx, field_b_ny)]) * metric_b[j];
        break;
    case advance_hy:
        curl = (field_a[index3(i + 1, j, k, field_a_nx, field_a_ny)]
                - field_a[index3(i, j, k, field_a_nx, field_a_ny)]) * metric_a[i]
            - (field_b[index3(i, j, k + 1, field_b_nx, field_b_ny)]
                - field_b[index3(i, j, k, field_b_nx, field_b_ny)]) * metric_b[k];
        break;
    default:
        curl = (field_a[index3(i, j + 1, k, field_a_nx, field_a_ny)]
                - field_a[index3(i, j, k, field_a_nx, field_a_ny)]) * metric_a[j]
            - (field_b[index3(i + 1, j, k, field_b_nx, field_b_ny)]
                - field_b[index3(i, j, k, field_b_nx, field_b_ny)]) * metric_b[i];
        break;
    }

    output[output_index] = coefficient_a[medium] * old_value + coefficient_b[medium] * curl;
}

template <typename Real>
__global__ void update_cpml_kernel(const CpmlJob* jobs, int job_count, int max_blocks) {
    const int job_index = blockIdx.y;
    if (job_index >= job_count) {
        return;
    }
    const CpmlJob& job = jobs[job_index];
    const int nx = job.region_xe - job.region_xi + 1;
    const int ny = job.region_ye - job.region_yi + 1;
    const int nz = job.region_ze - job.region_zi + 1;
    const int total = nx * ny * nz;
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= total || blockIdx.x >= max_blocks) {
        return;
    }

    const int global_i = job.region_xi + index % nx;
    const int global_j = job.region_yi + (index / nx) % ny;
    const int global_k = job.region_zi + index / (nx * ny);
    const int output_i = global_i - job.output.lo_x;
    const int output_j = global_j - job.output.lo_y;
    const int output_k = global_k - job.output.lo_z;
    const int media_i = global_i - job.media.lo_x;
    const int media_j = global_j - job.media.lo_y;
    const int media_k = global_k - job.media.lo_z;
    const int source_i = global_i - job.source.lo_x;
    const int source_j = global_j - job.source.lo_y;
    const int source_k = global_k - job.source.lo_z;
    const int psi_i = global_i - job.psi.lo_x;
    const int psi_j = global_j - job.psi.lo_y;
    const int psi_k = global_k - job.psi.lo_z;

    int source_neighbor_i = source_i;
    int source_neighbor_j = source_j;
    int source_neighbor_k = source_k;
    int coefficient_index;
    if (job.axis == 1) {
        source_neighbor_i += job.source_direction;
        coefficient_index = global_i;
    } else if (job.axis == 2) {
        source_neighbor_j += job.source_direction;
        coefficient_index = global_j;
    } else {
        source_neighbor_k += job.source_direction;
        coefficient_index = global_k;
    }

    const auto* source = static_cast<const Real*>(job.source.data);
    const auto* media = static_cast<const std::int16_t*>(job.media.data);
    auto* psi = static_cast<Real*>(job.psi.data);
    auto* output = static_cast<Real*>(job.output.data);
    const auto* coefficient_b = static_cast<const Real*>(job.coefficient_b.data);
    const auto* coefficient_c = static_cast<const Real*>(job.coefficient_c.data);
    const auto* coefficient_g = static_cast<const Real*>(job.coefficient_g.data);
    const auto psi_index = index3(psi_i, psi_j, psi_k, job.psi.nx, job.psi.ny);
    const auto output_index = index3(output_i, output_j, output_k, job.output.nx, job.output.ny);
    const auto medium_index = index3(media_i, media_j, media_k, job.media.nx, job.media.ny);
    const auto source_index = index3(source_i, source_j, source_k, job.source.nx, job.source.ny);
    const auto source_neighbor_index = index3(
            source_neighbor_i, source_neighbor_j, source_neighbor_k,
            job.source.nx, job.source.ny);
    const int medium = media[medium_index];
    const int coefficient_offset = coefficient_index - job.coefficient_b.lo;
    const Real source_difference = source[source_neighbor_index] - source[source_index];

    psi[psi_index] = coefficient_b[coefficient_offset] * psi[psi_index]
        + coefficient_c[coefficient_index - job.coefficient_c.lo] * source_difference;
    atomicAdd(output + output_index, static_cast<Real>(job.update_sign)
            * coefficient_g[medium - job.coefficient_g.lo] * psi[psi_index]);
}

template <typename Real>
__global__ void update_pml_body_kernel(
        int count,
        Real* psi,
        const Real* coefficient_b,
        const Real* coefficient_c,
        const Real* difference,
        const Real* coefficient_g,
        Real* field,
        int update_sign) {
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) {
        return;
    }
    psi[index] = coefficient_b[index] * psi[index] + coefficient_c[index] * difference[index];
    field[index] += static_cast<Real>(update_sign) * coefficient_g[index] * psi[index];
}

bool has_extent(const FieldView& view, int max_i, int max_j, int max_k) {
    return view.data != nullptr && view.nx > max_i && view.ny > max_j && view.nz > max_k;
}

struct DeviceAllocation {
    DeviceBuffer buffer;
    const void* host_data;
    std::size_t bytes;
    bool copy_back;
};

class DeviceAllocationMap {
public:
    cudaError_t get(const void* host_data, std::size_t bytes, bool copy_back, void*& device_data) {
        if (host_data == nullptr || bytes == 0) {
            return cudaErrorInvalidValue;
        }
        const auto found = allocations_.find(host_data);
        if (found != allocations_.end()) {
            auto& allocation = *found->second;
            if (allocation.bytes != bytes) {
                return cudaErrorInvalidValue;
            }
            allocation.copy_back = allocation.copy_back || copy_back;
            device_data = allocation.buffer.as<void>();
            return cudaSuccess;
        }

        auto allocation = std::make_unique<DeviceAllocation>();
        allocation->host_data = host_data;
        allocation->bytes = bytes;
        allocation->copy_back = copy_back;
        const auto status = allocation->buffer.copy_from(host_data, bytes);
        if (status != cudaSuccess) {
            return status;
        }
        device_data = allocation->buffer.as<void>();
        allocations_.emplace(host_data, allocation.get());
        allocations_storage_.push_back(std::move(allocation));
        return cudaSuccess;
    }

    cudaError_t copy_back() const {
        for (const auto& allocation : allocations_storage_) {
            if (!allocation->copy_back) {
                continue;
            }
            const auto status = allocation->buffer.copy_to(
                    const_cast<void*>(allocation->host_data), allocation->bytes);
            if (status != cudaSuccess) {
                return status;
            }
        }
        return cudaSuccess;
    }

private:
    std::unordered_map<const void*, DeviceAllocation*> allocations_;
    std::vector<std::unique_ptr<DeviceAllocation>> allocations_storage_;
};

bool field_bytes(const FieldView& view, std::size_t element_size, std::size_t& bytes) {
    if (view.data == nullptr || view.nx <= 0 || view.ny <= 0 || view.nz <= 0) {
        return false;
    }
    const auto max = std::numeric_limits<std::size_t>::max();
    auto elements = static_cast<std::size_t>(view.nx);
    if (elements > max / static_cast<std::size_t>(view.ny)) return false;
    elements *= static_cast<std::size_t>(view.ny);
    if (elements > max / static_cast<std::size_t>(view.nz)) return false;
    elements *= static_cast<std::size_t>(view.nz);
    if (elements > max / element_size) return false;
    bytes = elements * element_size;
    return true;
}

bool vector_bytes(const VectorView& view, std::size_t element_size, std::size_t& bytes) {
    if (view.data == nullptr || view.length <= 0) {
        return false;
    }
    const auto elements = static_cast<std::size_t>(view.length);
    if (elements > std::numeric_limits<std::size_t>::max() / element_size) {
        return false;
    }
    bytes = elements * element_size;
    return true;
}

template <typename Real>
cudaError_t prepare_cpml_job(
        const CpmlJob& job,
        DeviceAllocationMap& allocations,
        CpmlJob& device_job,
        int& total);

template <typename Real>
cudaError_t run_cpml_jobs(const CpmlJob* jobs, int job_count) {
    if (job_count == 0) return cudaSuccess;

    DeviceAllocationMap allocations;
    std::vector<CpmlJob> device_jobs(static_cast<std::size_t>(job_count));
    int max_blocks = 0;
    constexpr int threads = 256;
    for (int i = 0; i < job_count; ++i) {
        int total;
        const auto status = prepare_cpml_job<Real>(
                jobs[i], allocations, device_jobs[i], total);
        if (status != cudaSuccess) return status;
        max_blocks = std::max(max_blocks, (total - 1) / threads + 1);
    }

    DeviceBuffer device_job_buffer;
    const auto job_bytes = device_jobs.size() * sizeof(CpmlJob);
    auto status = device_job_buffer.copy_from(device_jobs.data(), job_bytes);
    if (status != cudaSuccess) return status;

    const dim3 grid(static_cast<unsigned int>(max_blocks),
            static_cast<unsigned int>(job_count));
    update_cpml_kernel<Real><<<grid, threads>>>(
            device_job_buffer.as<CpmlJob>(), job_count, max_blocks);
    status = cudaGetLastError();
    if (status != cudaSuccess) return status;
    status = cudaDeviceSynchronize();
    if (status != cudaSuccess) return status;
    return allocations.copy_back();
}

template <typename Real>
cudaError_t run_pml_body(int count, const VectorView& psi, const VectorView& coefficient_b,
        const VectorView& coefficient_c, const VectorView& difference,
        const VectorView& coefficient_g, const VectorView& field, int update_sign) {
    if (count <= 0) return cudaSuccess;
    if (update_sign != -1 && update_sign != 1) return cudaErrorInvalidValue;
    const VectorView* vectors[] = {
        &psi, &coefficient_b, &coefficient_c, &difference, &coefficient_g, &field
    };
    for (const auto* vector : vectors) {
        if (vector->data == nullptr || vector->length != count) {
            return cudaErrorInvalidValue;
        }
    }

    DeviceBuffer device_psi, device_coefficient_b, device_coefficient_c;
    DeviceBuffer device_difference, device_coefficient_g, device_field;
    const auto bytes = static_cast<std::size_t>(count) * sizeof(Real);
    auto status = device_psi.copy_from(psi.data, bytes);
    if (status != cudaSuccess) return status;
    status = device_coefficient_b.copy_from(coefficient_b.data, bytes);
    if (status != cudaSuccess) return status;
    status = device_coefficient_c.copy_from(coefficient_c.data, bytes);
    if (status != cudaSuccess) return status;
    status = device_difference.copy_from(difference.data, bytes);
    if (status != cudaSuccess) return status;
    status = device_coefficient_g.copy_from(coefficient_g.data, bytes);
    if (status != cudaSuccess) return status;
    status = device_field.copy_from(field.data, bytes);
    if (status != cudaSuccess) return status;

    constexpr int threads = 256;
    const int blocks = (count - 1) / threads + 1;
    update_pml_body_kernel<Real><<<blocks, threads>>>(
            count,
            device_psi.as<Real>(),
            device_coefficient_b.as<Real>(),
            device_coefficient_c.as<Real>(),
            device_difference.as<Real>(),
            device_coefficient_g.as<Real>(),
            device_field.as<Real>(),
            update_sign);
    status = cudaGetLastError();
    if (status != cudaSuccess) return status;
    status = cudaDeviceSynchronize();
    if (status != cudaSuccess) return status;
    status = device_psi.copy_to(psi.data, bytes);
    if (status != cudaSuccess) return status;
    return device_field.copy_to(field.data, bytes);
}

template <typename Real>
cudaError_t prepare_cpml_job(
        const CpmlJob& job,
        DeviceAllocationMap& allocations,
        CpmlJob& device_job,
        int& total) {
    const int nx = job.region_xe - job.region_xi + 1;
    const int ny = job.region_ye - job.region_yi + 1;
    const int nz = job.region_ze - job.region_zi + 1;
    if (nx <= 0 || ny <= 0 || nz <= 0
            || nx > (std::numeric_limits<int>::max() / ny) / nz
            || (job.axis < 1 || job.axis > 3)
            || (job.source_direction != -1 && job.source_direction != 1)
            || (job.update_sign != -1 && job.update_sign != 1)) {
        return cudaErrorInvalidValue;
    }
    total = nx * ny * nz;

    const int region_max_i = job.region_xe;
    const int region_max_j = job.region_ye;
    const int region_max_k = job.region_ze;
    int source_min_i = job.region_xi, source_max_i = region_max_i;
    int source_min_j = job.region_yi, source_max_j = region_max_j;
    int source_min_k = job.region_zi, source_max_k = region_max_k;
    if (job.axis == 1) {
        source_min_i += job.source_direction < 0 ? -1 : 0;
        source_max_i += job.source_direction > 0 ? 1 : 0;
    } else if (job.axis == 2) {
        source_min_j += job.source_direction < 0 ? -1 : 0;
        source_max_j += job.source_direction > 0 ? 1 : 0;
    } else {
        source_min_k += job.source_direction < 0 ? -1 : 0;
        source_max_k += job.source_direction > 0 ? 1 : 0;
    }

    if (!has_extent(job.output, region_max_i - job.output.lo_x,
                region_max_j - job.output.lo_y, region_max_k - job.output.lo_z)
            || job.region_xi < job.output.lo_x || job.region_yi < job.output.lo_y
            || job.region_zi < job.output.lo_z
            || region_max_i < job.output.lo_x || region_max_j < job.output.lo_y
            || region_max_k < job.output.lo_z
            || !has_extent(job.media, region_max_i - job.media.lo_x,
                region_max_j - job.media.lo_y, region_max_k - job.media.lo_z)
            || job.region_xi < job.media.lo_x || job.region_yi < job.media.lo_y
            || job.region_zi < job.media.lo_z
            || region_max_i < job.media.lo_x || region_max_j < job.media.lo_y
            || region_max_k < job.media.lo_z
            || !has_extent(job.psi, region_max_i - job.psi.lo_x,
                region_max_j - job.psi.lo_y, region_max_k - job.psi.lo_z)
            || job.region_xi < job.psi.lo_x || job.region_yi < job.psi.lo_y
            || job.region_zi < job.psi.lo_z
            || region_max_i < job.psi.lo_x || region_max_j < job.psi.lo_y
            || region_max_k < job.psi.lo_z
            || source_min_i < job.source.lo_x || source_min_j < job.source.lo_y
            || source_min_k < job.source.lo_z
            || source_max_i >= job.source.lo_x + job.source.nx
            || source_max_j >= job.source.lo_y + job.source.ny
            || source_max_k >= job.source.lo_z + job.source.nz
            || job.coefficient_g.lo > 0
            || job.coefficient_g.length <= 0) {
        return cudaErrorInvalidValue;
    }

    const int coefficient_min = job.axis == 1 ? job.region_xi
        : job.axis == 2 ? job.region_yi : job.region_zi;
    const int coefficient_max = job.axis == 1 ? job.region_xe
        : job.axis == 2 ? job.region_ye : job.region_ze;
    if (coefficient_min < job.coefficient_b.lo || coefficient_min < job.coefficient_c.lo
            || coefficient_max >= job.coefficient_b.lo + job.coefficient_b.length
            || coefficient_max >= job.coefficient_c.lo + job.coefficient_c.length) {
        return cudaErrorInvalidValue;
    }

    std::size_t bytes;
    void* device_data;
    if (!field_bytes(job.output, sizeof(Real), bytes)) return cudaErrorInvalidValue;
    auto status = allocations.get(job.output.data, bytes, true, device_data);
    if (status != cudaSuccess) return status;
    device_job.output = job.output;
    device_job.output.data = device_data;

    if (!field_bytes(job.media, sizeof(std::int16_t), bytes)) return cudaErrorInvalidValue;
    status = allocations.get(job.media.data, bytes, false, device_data);
    if (status != cudaSuccess) return status;
    device_job.media = job.media;
    device_job.media.data = device_data;

    if (!field_bytes(job.source, sizeof(Real), bytes)) return cudaErrorInvalidValue;
    status = allocations.get(job.source.data, bytes, false, device_data);
    if (status != cudaSuccess) return status;
    device_job.source = job.source;
    device_job.source.data = device_data;

    if (!field_bytes(job.psi, sizeof(Real), bytes)) return cudaErrorInvalidValue;
    status = allocations.get(job.psi.data, bytes, true, device_data);
    if (status != cudaSuccess) return status;
    device_job.psi = job.psi;
    device_job.psi.data = device_data;

    const VectorView* host_vectors[] = {
        &job.coefficient_b, &job.coefficient_c, &job.coefficient_g
    };
    VectorView* device_vectors[] = {
        &device_job.coefficient_b, &device_job.coefficient_c, &device_job.coefficient_g
    };
    for (int i = 0; i < 3; ++i) {
        if (!vector_bytes(*host_vectors[i], sizeof(Real), bytes)) {
            return cudaErrorInvalidValue;
        }
        status = allocations.get(host_vectors[i]->data, bytes, false, device_data);
        if (status != cudaSuccess) return status;
        *device_vectors[i] = *host_vectors[i];
        device_vectors[i]->data = device_data;
    }
    device_job.region_xi = job.region_xi;
    device_job.region_xe = job.region_xe;
    device_job.region_yi = job.region_yi;
    device_job.region_ye = job.region_ye;
    device_job.region_zi = job.region_zi;
    device_job.region_ze = job.region_ze;
    device_job.axis = job.axis;
    device_job.source_direction = job.source_direction;
    device_job.update_sign = job.update_sign;
    return cudaSuccess;
}

bool has_extent(const VectorView& view, int max_index) {
    return view.data != nullptr && view.length > max_index;
}

template <typename Real>
cudaError_t copy_real_field(DeviceBuffer& device, const FieldView& host) {
    const auto elements = static_cast<std::size_t>(host.nx) * host.ny * host.nz;
    if (elements > std::numeric_limits<std::size_t>::max() / sizeof(Real)) {
        return cudaErrorInvalidValue;
    }
    return device.copy_from(host.data, elements * sizeof(Real));
}

template <typename Real>
cudaError_t copy_real_vector(DeviceBuffer& device, const VectorView& host) {
    const auto elements = static_cast<std::size_t>(host.length);
    if (elements > std::numeric_limits<std::size_t>::max() / sizeof(Real)) {
        return cudaErrorInvalidValue;
    }
    return device.copy_from(host.data, elements * sizeof(Real));
}

template <typename Real>
cudaError_t run_update(
        int operation,
        const FieldView& output,
        const FieldView& media,
        const FieldView& field_a,
        const FieldView& field_b,
        const VectorView& metric_a,
        const VectorView& metric_b,
        const VectorView& coefficient_a,
        const VectorView& coefficient_b,
        const Sweep& sweep) {
    if (sweep.nx <= 0 || sweep.ny <= 0 || sweep.nz <= 0) {
        return cudaSuccess;
    }
    if (sweep.nx > (std::numeric_limits<int>::max() / sweep.ny) / sweep.nz) {
        return cudaErrorInvalidValue;
    }

    int max_a_i = sweep.nx, max_a_j = sweep.ny, max_a_k = sweep.nz;
    int max_b_i = sweep.nx, max_b_j = sweep.ny, max_b_k = sweep.nz;
    int metric_a_index = 0, metric_b_index = 0;
    switch (operation) {
    case advance_ex:
        metric_a_index = sweep.ny;
        metric_b_index = sweep.nz;
        break;
    case advance_ey:
        metric_a_index = sweep.nz;
        metric_b_index = sweep.nx;
        break;
    case advance_ez:
        metric_a_index = sweep.nx;
        metric_b_index = sweep.ny;
        break;
    case advance_hx:
        max_a_k++;
        max_b_j++;
        metric_a_index = sweep.nz;
        metric_b_index = sweep.ny;
        break;
    case advance_hy:
        max_a_i++;
        max_b_k++;
        metric_a_index = sweep.nx;
        metric_b_index = sweep.nz;
        break;
    case advance_hz:
        max_a_j++;
        max_b_i++;
        metric_a_index = sweep.ny;
        metric_b_index = sweep.nx;
        break;
    default:
        return cudaErrorInvalidValue;
    }

    if (!has_extent(output, sweep.nx, sweep.ny, sweep.nz)
            || !has_extent(media, sweep.nx, sweep.ny, sweep.nz)
            || !has_extent(field_a, max_a_i, max_a_j, max_a_k)
            || !has_extent(field_b, max_b_i, max_b_j, max_b_k)
            || !has_extent(metric_a, metric_a_index)
            || !has_extent(metric_b, metric_b_index)
            || coefficient_a.data == nullptr || coefficient_b.data == nullptr
            || coefficient_a.length <= 0 || coefficient_b.length <= 0) {
        return cudaErrorInvalidValue;
    }

    DeviceBuffer device_output, device_media, device_a, device_b;
    DeviceBuffer device_metric_a, device_metric_b, device_coefficient_a, device_coefficient_b;
    cudaError_t status;

#define COPY_REAL_FIELD(device, host) \
    do { status = copy_real_field<Real>((device), (host)); if (status != cudaSuccess) return status; } while (false)
#define COPY_REAL_VECTOR(device, host) \
    do { status = copy_real_vector<Real>((device), (host)); if (status != cudaSuccess) return status; } while (false)

    COPY_REAL_FIELD(device_output, output);
    status = device_media.copy_from(media.data,
            static_cast<std::size_t>(media.nx) * media.ny * media.nz * sizeof(std::int16_t));
    if (status != cudaSuccess) return status;
    COPY_REAL_FIELD(device_a, field_a);
    COPY_REAL_FIELD(device_b, field_b);
    COPY_REAL_VECTOR(device_metric_a, metric_a);
    COPY_REAL_VECTOR(device_metric_b, metric_b);
    COPY_REAL_VECTOR(device_coefficient_a, coefficient_a);
    COPY_REAL_VECTOR(device_coefficient_b, coefficient_b);

    const int total = sweep.nx * sweep.ny * sweep.nz;
    constexpr int threads = 256;
    const int blocks = (total - 1) / threads + 1;
    update_field_kernel<Real><<<blocks, threads>>>(
            operation,
            device_output.as<Real>(),
            device_media.as<std::int16_t>(),
            device_a.as<Real>(),
            device_b.as<Real>(),
            device_metric_a.as<Real>(),
            device_metric_b.as<Real>(),
            device_coefficient_a.as<Real>(),
            device_coefficient_b.as<Real>(),
            output.nx,
            output.ny,
            media.nx,
            media.ny,
            field_a.nx,
            field_a.ny,
            field_b.nx,
            field_b.ny,
            sweep.nx,
            sweep.ny,
            sweep.nz,
            total);
    status = cudaGetLastError();
    if (status != cudaSuccess) return status;
    status = cudaDeviceSynchronize();
    if (status != cudaSuccess) return status;
    status = device_output.copy_to(output.data,
            static_cast<std::size_t>(output.nx) * output.ny * output.nz * sizeof(Real));
    if (status != cudaSuccess) return status;

#undef COPY_REAL_FIELD
#undef COPY_REAL_VECTOR
    return cudaSuccess;
}

cudaError_t report_cuda_status(cudaError_t status, const char* operation) {
    if (status != cudaSuccess) {
        std::fprintf(stderr, "%s failed: %s\n", operation, cudaGetErrorString(status));
    }
    return status;
}

} // namespace

extern "C" int semba_cuda_advance_cpml_jobs(const CpmlJob* jobs, int job_count) {
    if (jobs == nullptr || job_count < 0) {
        return static_cast<int>(report_cuda_status(
                cudaErrorInvalidValue, "CUDA CPML update"));
    }
#ifdef SEMBA_CUDA_DOUBLE_PRECISION
    return static_cast<int>(report_cuda_status(
            run_cpml_jobs<double>(jobs, job_count), "CUDA CPML update"));
#else
    return static_cast<int>(report_cuda_status(
            run_cpml_jobs<float>(jobs, job_count), "CUDA CPML update"));
#endif
}

extern "C" int semba_cuda_advance_pml_body(
        int count,
        const VectorView* psi,
        const VectorView* coefficient_b,
        const VectorView* coefficient_c,
        const VectorView* difference,
        const VectorView* coefficient_g,
        const VectorView* field,
        int update_sign) {
    if (psi == nullptr || coefficient_b == nullptr || coefficient_c == nullptr
            || difference == nullptr || coefficient_g == nullptr || field == nullptr) {
        return static_cast<int>(report_cuda_status(
                cudaErrorInvalidValue, "CUDA PML-body update"));
    }
#ifdef SEMBA_CUDA_DOUBLE_PRECISION
    return static_cast<int>(report_cuda_status(
            run_pml_body<double>(count, *psi, *coefficient_b, *coefficient_c,
                    *difference, *coefficient_g, *field, update_sign),
            "CUDA PML-body update"));
#else
    return static_cast<int>(report_cuda_status(
            run_pml_body<float>(count, *psi, *coefficient_b, *coefficient_c,
                    *difference, *coefficient_g, *field, update_sign),
            "CUDA PML-body update"));
#endif
}

extern "C" int semba_cuda_advance_field(
        int operation,
        const FieldView* output,
        const FieldView* media,
        const FieldView* field_a,
        const FieldView* field_b,
        const VectorView* metric_a,
        const VectorView* metric_b,
        const VectorView* coefficient_a,
        const VectorView* coefficient_b,
        const Sweep* sweep) {
    if (output == nullptr || media == nullptr || field_a == nullptr || field_b == nullptr
            || metric_a == nullptr || metric_b == nullptr || coefficient_a == nullptr
            || coefficient_b == nullptr || sweep == nullptr) {
        std::fprintf(stderr, "CUDA Yee update received a null argument.\n");
        return static_cast<int>(cudaErrorInvalidValue);
    }

#ifdef SEMBA_CUDA_DOUBLE_PRECISION
    const auto status = run_update<double>(
            operation, *output, *media, *field_a, *field_b, *metric_a, *metric_b,
            *coefficient_a, *coefficient_b, *sweep);
#else
    const auto status = run_update<float>(
            operation, *output, *media, *field_a, *field_b, *metric_a, *metric_b,
            *coefficient_a, *coefficient_b, *sweep);
#endif

    if (status != cudaSuccess) {
        std::fprintf(stderr, "CUDA Yee update failed: %s\n", cudaGetErrorString(status));
    }
    return static_cast<int>(status);
}
