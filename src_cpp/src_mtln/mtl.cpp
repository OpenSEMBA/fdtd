#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <optional>
#include <memory>

#include "mtln_types.h"
#include "mtl_m.h"
#include "multipolar_expansion_m.h"

namespace utils_m {
    std::vector<std::vector<double>> inv(const std::vector<std::vector<double>>& A);
    std::vector<std::vector<double>> element_wise_invert(int n, const std::vector<std::vector<double>>& x);
    std::vector<std::vector<double>> eye(int dim);
    std::vector<double> getEigenValues(const std::vector<std::vector<double>>& matrix);
}

using mtln_types_m::multipolar_expansion_t;
using mtln_types_m::segment_t;
using mtln_types_m::transfer_impedance_per_meter_t;

struct generator_t {};

// Constants from FDETYPES_m
// extern const double pi;
const double pi = 3.14159265358979323846;
// extern const double mu_vacuum;
// extern const double c_vacuum;
const double mu_vacuum = 4.0 * 3.14159265358979323846e-7;
const double c_vacuum = 299792458.0;
// extern const int RKIND_wires; // Assuming this maps to a specific float/double precision indicator or is just a tag
// extern const int RKIND;
// extern const int RKIND_TIEMPO;

// MPI related
#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
extern int REALSIZE;
extern int INTEGERSIZE;
#endif

namespace mtl_m {

using mtln_types_m::multipolar_expansion_t;
using mtln_types_m::segment_t;
using mtln_types_m::transfer_impedance_per_meter_t;

    // bool isSegmentZOriented(int j);
    // bool isSegmentNextToLayerEnd(int j, int z_end);
    // bool isSegmentBeforeLayerEnd(int j, int z_end);
    // bool isSegmentZPositive(int j);
    // bool isSegmentAfterLayerEnd(int j, int z_end);
    // bool isSegmentNextToLayerInit(int j, int z_init);
    // bool isSegmentBeforeLayerInit(int j, int z_init);
    // bool isSegmentAfterLayerInit(int j, int z_init);

    // Constants for MPI
#ifdef CompileWithMPI
    constexpr int COMM_SEND = 1;
    constexpr int COMM_RECV = -1;
    constexpr int COMM_NONE = 0;
    constexpr int COMM_FIELD = 1;
    constexpr int COMM_V = 2;
    constexpr int COMM_BOTH = 3;

#endif

    // Implementation of mtl_t methods

    void mtl_t::initLC(const std::vector<std::vector<double>>& lpul, const std::vector<std::vector<double>>& cpul) {
        int n = static_cast<int>(lpul.size());
        for (int i = 0; i < static_cast<int>(this->lpul.size()); ++i) {
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    this->lpul[i][r][c] = lpul[r][c];
                }
            }
        }
        for (int i = 0; i < static_cast<int>(this->cpul.size()); ++i) {
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    this->cpul[i][r][c] = cpul[r][c];
                }
            }
        }
    }

    void mtl_t::allocatePULMatrices() {
        int n = this->number_of_conductors;
        int sz_step = this->step_size.size();
        
        // Allocate with source 0.0
        this->lpul.assign(sz_step, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
        this->cpul.assign(sz_step + 1, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
        this->rpul.assign(sz_step, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
        this->gpul.assign(sz_step + 1, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
    }

    void mtl_t::computeLCParameters(const multipolar_expansion_t& multipolar_expansion) {
        int n = this->cpul[0].size();
        std::vector<std::vector<double>> ppul(n, std::vector<double>(n));
        
        for (int i = 0; i < static_cast<int>(this->segments.size()); ++i) {
            // Assuming getCellInductanceOnBox and getCellCapacitanceOnBox return 2D vectors
            // Note: The Fortran code assigns the result of these functions directly to 2D slices of 3D arrays
            // This implies the functions return a 2D matrix.
            this->lpul[i] = multipolar_expansion_m::getCellInductanceOnBox(multipolar_expansion, this->segments[i].dualBox);
            this->cpul[i] = multipolar_expansion_m::getCellCapacitanceOnBox(multipolar_expansion, this->segments[i].dualBox);
            
            ppul = utils_m::element_wise_invert(n, this->cpul[i]);
            this->cpul[i] = utils_m::inv(ppul);
        }
        
        // Copy last slice to the extra slice
        int last_idx = static_cast<int>(this->segments.size()) - 1;
        this->cpul[this->segments.size()] = this->cpul[last_idx];
    }

    void mtl_t::computeLCParametersFromRadius(double rad) {
        double invMu = 1.0 / mu_vacuum;
        for (int i = 0; i < static_cast<int>(this->segments.size()); ++i) {
            double d1 = this->segments[i].d1;
            double d2 = this->segments[i].d2;
            
            // Fortran: (1.0_RKIND_wires / (4.0_RKIND_wires * pi*invMu))*(...)
            // Assuming RKIND_wires is just a kind specifier for double in C++ context here
            double term1 = (1.0 / (4.0 * pi * invMu)) * 
                (std::log((d1*d1 + d2*d2) / (4.0 * rad*rad)) + 
                 (d1/d2)*std::atan(d2/d1) + 
                 (d2/d1)*std::atan(d1/d2) + 
                 (pi * rad*rad) / (d2*d1) - 3.0);
            
            this->lpul[i] = std::vector<std::vector<double>>(this->number_of_conductors, 
                                                              std::vector<double>(this->number_of_conductors, term1));

            if ((rad < 0.3 * d1) || (rad < 0.3 * d2)) {
                double correction = 0.57 / (4.0 * pi * invMu);
                for(int r=0; r<this->number_of_conductors; ++r)
                    for(int c=0; c<this->number_of_conductors; ++c)
                        this->lpul[i][r][c] -= correction;
            }

            if ((rad > 0.3 * d1) || (rad > 0.3 * d2)) {
                double divisor = 1.0 - (pi * rad*rad) / (d1*d2);
                for(int r=0; r<this->number_of_conductors; ++r)
                    for(int c=0; c<this->number_of_conductors; ++c)
                        this->lpul[i][r][c] /= divisor;
            }
            
            // cpul = 1 / (lpul * c_vacuum^2)
            double c_v_sq = c_vacuum * c_vacuum;
            for(int r=0; r<this->number_of_conductors; ++r)
                for(int c=0; c<this->number_of_conductors; ++c)
                    this->cpul[i][r][c] = 1.0 / (this->lpul[i][r][c] * c_v_sq);
        }
        
        int last_idx = static_cast<int>(this->segments.size()) - 1;
        this->cpul[this->segments.size()] = this->cpul[last_idx];
    }

    void mtl_t::initDirections() {
        int n = this->number_of_conductors;
        int sz = this->step_size.size();
        
        // Fortran: reshape(source = [(this%step_size(j)*eye(this%number_of_conductors) , j = 1, size(this%step_size))], ... order=[2,3,1])
        // This creates a 3D array where the first dimension is step_size, second is conductor, third is conductor.
        // The source is a list of matrices. The reshape with order=[2,3,1] means the memory layout is filled such that
        // the first index varies slowest? No, Fortran order=[2,3,1] means the first index of the result corresponds to the 2nd index of the source list?
        // Actually, reshape fills the array in column-major order. The source is a list of matrices.
        // Let's interpret: res(i,j,k) comes from source(k) ?
        // If order=[2,3,1], then res(1,1,1) is source(1)(1,1). res(2,1,1) is source(1)(1,1) again? No.
        // Standard reshape: elements are taken from source in order and placed into result in order.
        // If source is a list of matrices M_1, M_2, ..., M_N.
        // And result shape is [N, n, n].
        // Then res(1,1,1) = M_1(1,1). res(1,2,1) = M_1(2,1). ... res(1,n,n) = M_1(n,n).
        // res(2,1,1) = M_2(1,1).
        // This matches the loop structure usually implied by such assignments in Fortran if not using explicit reshape.
        // However, the code uses reshape.
        // Let's just assign directly to match the logic:
        // this%du(j, :, :) = step_size(j) * eye(n)
        // But the shape is [size(step_size), n, n].
        // So du(i, r, c) = step_size(i) * (r==c ? 1 : 0)
        
        this->du.assign(sz, std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0)));
        for (int j = 0; j < sz; ++j) {
            double val = this->step_size[j];
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    if (r == c) {
                        this->du[j][r][c] = val;
                    } else {
                        this->du[j][r][c] = 0.0;
                    }
                }
            }
        }
    }

    void mtl_t::checkTimeStep(bool getMax, std::optional<double> dt) {
        if (dt.has_value()) {
            if (getMax) {
                double max_dt = this->getMaxTimeStep();
                if (dt.value() > max_dt && max_dt > 0) {
                    this->dt = max_dt;
                    std::cout << "dt larger than maximum permitted. Changed to dt = " << max_dt << std::endl;
                } else {
                    this->dt = dt.value();
                }
            } else {
                this->dt = dt.value();
            }
        } else {
            if (getMax) {
                this->dt = this->getMaxTimeStep();
            }
        }
    }

    void mtl_t::initRG(const std::vector<std::vector<double>>& rpul, const std::vector<std::vector<double>>& gpul) {
        int n = static_cast<int>(rpul.size());
        for (int i = 0; i < static_cast<int>(this->rpul.size()); ++i) {
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    this->rpul[i][r][c] = rpul[r][c];
                }
            }
        }
        for (int i = 0; i < static_cast<int>(this->gpul.size()); ++i) {
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    this->gpul[i][r][c] = gpul[r][c];
                }
            }
        }
    }

    double mtl_t::getMaxTimeStep() {
        // minval(pack(this%du, this%du /= 0))
        // This finds the minimum non-zero value in du
        double min_val = 1e30;
        bool found = false;
        for (const auto& layer : this->du) {
            for (const auto& row : layer) {
                for (double val : row) {
                    if (val != 0.0) {
                        if (val < min_val) {
                            min_val = val;
                            found = true;
                        }
                    }
                }
            }
        }
        
        if (!found) return 0.0; // Or handle error

        // maxval(this%getPhaseVelocities())
        auto phase_vels = this->getPhaseVelocities();
        double max_vel = 0.0;
        for (const auto& row : phase_vels) {
            for (double v : row) {
                if (v > max_vel) max_vel = v;
            }
        }
        
        if (max_vel == 0.0) return 0.0;

        return min_val / max_vel;
    }

    std::vector<std::vector<double>> mtl_t::getPhaseVelocities() {
        int n = this->number_of_conductors;
        int sz = this->step_size.size();
        std::vector<std::vector<double>> res(sz, std::vector<double>(n));
        
        for (int k = 0; k < sz; ++k) {
            // matmul(lpul(k), cpul(k+1))
            // Note: cpul has size sz+1, so k+1 is valid for k < sz
            std::vector<std::vector<double>> matmul_res(n, std::vector<double>(n, 0.0));
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    double sum = 0.0;
                    for (int l = 0; l < n; ++l) {
                        sum += this->lpul[k][i][l] * this->cpul[k+1][l][j];
                    }
                    matmul_res[i][j] = sum;
                }
            }
            
            // getEigenValues returns a vector. The code uses ev(1:this%number_of_conductors)
            // Assuming getEigenValues returns eigenvalues of the matrix.
            // The matrix is N x N. It should return N eigenvalues.
            // The Fortran code calls getEigenValues(dble(...)).
            std::vector<double> ev = utils_m::getEigenValues(matmul_res);
            
            // res(k,:) = 1.0/sqrt(ev(1:n))
            for (int i = 0; i < n; ++i) {
                if (i < static_cast<int>(ev.size())) {
                    res[k][i] = 1.0 / std::sqrt(ev[i]);
                } else {
                    res[k][i] = 0.0; // Fallback
                }
            }
        }
        return res;
    }

    void mtl_t::setTimeStep(int numberOfSteps, double finalTime) {
        this->dt = finalTime / numberOfSteps;
    }

#ifdef CompileWithMPI
    void mtl_t::initStepSizeAndFieldSegments(const std::vector<double>& step_size, const std::vector<segment_t>& segments, const std::vector<std::vector<int>>& layer_indices) {
        int n = 0;
        for (int j = 0; j < static_cast<int>(layer_indices.size()); ++j) {
            n += layer_indices[j][1] - layer_indices[j][0] + 1;
        }
        n += static_cast<int>(layer_indices.size()) - 1;

        this->step_size.resize(n);
        this->segments.resize(n);

        int idx = 0;
        for (int j = 0; j < static_cast<int>(layer_indices.size()); ++j) {
            int start = layer_indices[j][0];
            int end = layer_indices[j][1];
            int count = end - start + 1;
            
            // Copy step_size
            for (int k = 0; k < count; ++k) {
                this->step_size[idx + k] = step_size[start + k];
            }
            // Copy segments
            for (int k = 0; k < count; ++k) {
                this->segments[idx + k] = segments[start + k];
            }

            if (j != static_cast<int>(layer_indices.size()) - 1) {
                // Duplicate the last element and set orientation to -1
                this->step_size[idx + count] = this->step_size[idx + count - 1];
                this->segments[idx + count] = this->segments[idx + count - 1];
                this->segments[idx + count].orientation = -1;
                idx += count + 1;
            } else {
                idx += count;
            }
        }
    }

    void mtl_t::initCommunicators(const std::vector<int>& alloc_z) {
        int rank, ierr;
        MPI_Comm_rank(SUBCOMM_MPI, &rank, &ierr);
        this->mpi_comm.rank = rank;
        this->mpi_comm.comms.clear();
        
        int z_init = alloc_z[0];
        int z_end = alloc_z[1];

        for (int j = 0; j < static_cast<int>(this->segments.size()); ++j) {
            if (this->segments[j].orientation == -1) continue;

            int z = this->segments[j].z;
            
            // Check isSegmentZOriented(j)
            // Assuming this function is available in the namespace or global scope
            if (!isSegmentZOriented(j) && ((z == z_end) || (z == z_init + 1))) {
                int n = static_cast<int>(this->mpi_comm.comms.size());
                std::vector<communicator_t> aux_comm = this->mpi_comm.comms;
                aux_comm.resize(n + 1);
                
                aux_comm[n].field_index = j;
                aux_comm[n].comm_type = COMM_FIELD;
                aux_comm[n].v_index = -1;
                
                if (z == z_end) {
                    aux_comm[n].delta_rank = 1;
                    aux_comm[n].comm_task = COMM_RECV;
                } else if (z == z_init + 1) {
                    aux_comm[n].delta_rank = -1;
                    aux_comm[n].comm_task = COMM_SEND;
                }
                
                this->mpi_comm.comms = aux_comm;
            }

            // Check isSegmentZOriented(j) AND (isSegmentNextToLayerEnd OR isSegmentNextToLayerInit)
            if (isSegmentZOriented(j) && (isSegmentNextToLayerEnd(j, z_end) || isSegmentNextToLayerInit(j, z_init))) {
                int n = static_cast<int>(this->mpi_comm.comms.size());
                std::vector<communicator_t> aux_comm = this->mpi_comm.comms;
                aux_comm.resize(n + 1);
                
                aux_comm[n].field_index = j;
                aux_comm[n].comm_type = COMM_BOTH;

                if (isSegmentNextToLayerEnd(j, z_end)) {
                    aux_comm[n].delta_rank = 1;
                    if (isSegmentBeforeLayerEnd(j, z_end)) {
                        aux_comm[n].comm_task = COMM_SEND;
                        if (isSegmentZPositive(j)) aux_comm[n].v_index = j;
                        else aux_comm[n].v_index = j + 1;
                    } else if (isSegmentAfterLayerEnd(j, z_end)) {
                        aux_comm[n].comm_task = COMM_RECV;
                        if (isSegmentZPositive(j)) aux_comm[n].v_index = j + 1;
                        else aux_comm[n].v_index = j;
                    }
                } else if (isSegmentNextToLayerInit(j, z_init)) {
                    aux_comm[n].delta_rank = -1;
                    if (isSegmentBeforeLayerInit(j, z_init)) {
                        aux_comm[n].comm_task = COMM_RECV;
                        if (isSegmentZPositive(j)) aux_comm[n].v_index = j;
                        else aux_comm[n].v_index = j + 1;
                    } else if (isSegmentAfterLayerInit(j, z_init)) {
                        aux_comm[n].comm_task = COMM_SEND;
                        if (isSegmentZPositive(j)) aux_comm[n].v_index = j + 1;
                        else aux_comm[n].v_index = j;
                    }
                }
                
                this->mpi_comm.comms = aux_comm;
            }
        }
    }
#endif

    // Interface function implementations

    mtl_t mtl_shielded(
        const std::vector<std::vector<double>>& lpul,
        const std::vector<std::vector<double>>& cpul,
        const std::vector<std::vector<double>>& rpul,
        const std::vector<std::vector<double>>& gpul,
        const std::vector<double>& step_size,
        const std::string& name,
        const std::vector<segment_t>& segments,
        double dt,
        const std::string& parent_name,
        int conductor_in_parent,
        const transfer_impedance_per_meter_t& transfer_impedance,
        std::optional<std::vector<std::vector<int>>> layer_indices,
        std::optional<bool> bundle_in_layer,
        std::optional<std::vector<int>> alloc_z
    ) {
        mtl_t res;
        
        checkPULDimensions(lpul, cpul, rpul, gpul);
        
        res.name = name;
        res.number_of_conductors = lpul.size();
        
#ifdef CompileWithMPI
        if (layer_indices.has_value()) {
            if (!layer_indices.value().empty()) {
                res.initStepSizeAndFieldSegments(step_size, segments, layer_indices.value());
                res.initCommunicators(alloc_z.value());
            } else {
                res.step_size.resize(0);
                res.segments.resize(0);
                res.mpi_comm.comms.resize(0);
            }
        } else {
            res.step_size = step_size;
            res.layer_indices.resize(0, 0);
            res.mpi_comm.comms.resize(0);
            res.segments = segments;
        }
#else
        res.step_size = step_size;
        res.segments = segments;
#endif

        res.initDirections();
        res.allocatePULMatrices();
        res.initLC(lpul, cpul);
        res.initRG(rpul, gpul);
        
        // checkTimeStep(getMax = (lpul(1,1) /= 0.0), dt = dt)
        // lpul(1,1) is lpul[0][0] in 0-based indexing
        bool getMax = (lpul[0][0] != 0.0);
        res.checkTimeStep(getMax, dt);
        
        res.parent_name = parent_name;
        res.conductor_in_parent = conductor_in_parent;
        res.transfer_impedance = transfer_impedance;
        res.lumped_elements = dispersive_m::lumped_t(res.number_of_conductors, 0,
                                                     static_cast<int>(res.step_size.size()), res.dt);
        
        return res;
    }

    mtl_t mtl_unshielded(
        const std::vector<std::vector<double>>& lpul,
        const std::vector<std::vector<double>>& cpul,
        const std::vector<std::vector<double>>& rpul,
        const std::vector<std::vector<double>>& gpul,
        const std::vector<double>& step_size,
        const std::string& name,
        const std::vector<segment_t>& segments,
        double dt,
        const std::vector<multipolar_expansion_t>& multipolar_expansion,
        double radius,
        std::optional<std::vector<std::vector<int>>> layer_indices,
        std::optional<bool> bundle_in_layer,
        std::optional<std::vector<int>> alloc_z
    ) {
        mtl_t res;
        
        checkPULDimensions(lpul, cpul, rpul, gpul);
        
        res.name = name;
        res.number_of_conductors = lpul.size();
        
#ifdef CompileWithMPI
        if (layer_indices.has_value()) {
            if (!layer_indices.value().empty()) {
                res.initStepSizeAndFieldSegments(step_size, segments, layer_indices.value());
                res.initCommunicators(alloc_z.value());
                res.layer_indices = layer_indices.value();
                res.bundle_in_layer = bundle_in_layer.value();
            } else {
                res.step_size.resize(0);
                res.segments.resize(0);
                res.mpi_comm.comms.resize(0);
            }
        } else {
            res.step_size = step_size;
            res.layer_indices.resize(0, 0);
            res.mpi_comm.comms.resize(0);
            res.segments = segments;
        }
#else
        res.step_size = step_size;
        res.segments = segments;
#endif

        res.initDirections();
        res.allocatePULMatrices();
        
        if (!res.step_size.empty()) {
            if (!multipolar_expansion.empty()) {
                res.computeLCParameters(multipolar_expansion[0]);
            } else if (radius != 0.0) {
                res.computeLCParametersFromRadius(radius);
            } else {
                res.initLC(lpul, cpul);
            }
        }
        
        res.initRG(rpul, gpul);
        
        bool getMax = (lpul[0][0] != 0.0);
        res.checkTimeStep(getMax, dt);
        
        res.lumped_elements = dispersive_m::lumped_t(res.number_of_conductors, 0,
                                                     static_cast<int>(res.step_size.size()), res.dt);
        
        return res;
    }

    void checkPULDimensions(
        const std::vector<std::vector<double>>& lpul,
        const std::vector<std::vector<double>>& cpul,
        const std::vector<std::vector<double>>& rpul,
        const std::vector<std::vector<double>>& gpul
    ) {
        int n_lpul = lpul.size();
        if (n_lpul != lpul[0].size()) {
            throw std::runtime_error("PUL matrices are not square");
        }
        
        int n_cpul = cpul.size();
        if (n_cpul != cpul[0].size()) {
            throw std::runtime_error("PUL matrices are not square");
        }
        
        int n_rpul = rpul.size();
        if (n_rpul != rpul[0].size()) {
            throw std::runtime_error("PUL matrices are not square");
        }
        
        int n_gpul = gpul.size();
        if (n_gpul != gpul[0].size()) {
            throw std::runtime_error("PUL matrices are not square");
        }
        
        if (n_lpul != n_cpul || n_lpul != n_rpul || n_lpul != n_gpul) {
            throw std::runtime_error("PUL matrices do not have the same dimensions");
        }
    }

} // namespace mtl_m