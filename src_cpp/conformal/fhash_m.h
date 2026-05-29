#ifndef FHASH_M_H
#define FHASH_M_H

#include <any>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

namespace fhash_m {

constexpr int FHASH_KEY_NOT_FOUND = -1;
constexpr int FHASH_EMPTY_TABLE = -3;
constexpr int FHASH_DEFAULT_ALLOCATION = 127;

constexpr std::uint64_t FNV_OFFSET_32 = 2166136261ULL;
constexpr std::uint64_t FNV_PRIME_32 = 16777619ULL;

inline std::uint64_t fnv_1a_int32(std::uint64_t seed, std::int32_t input) {
    std::uint64_t hash = seed;
    unsigned char bytes[4];
    std::memcpy(bytes, &input, sizeof(input));
    for (unsigned char byte : bytes) {
        hash ^= static_cast<std::uint64_t>(byte);
        hash *= FNV_PRIME_32;
    }
    return hash;
}

inline std::uint64_t fnv_1a_string(std::uint64_t seed, const std::string& input) {
    std::uint64_t hash = seed;
    for (unsigned char ch : input) {
        hash ^= static_cast<std::uint64_t>(ch);
        hash *= FNV_PRIME_32;
    }
    return hash;
}

inline std::uint64_t fnv_1a_int32_1d(const std::vector<std::int32_t>& input) {
    std::uint64_t hash = FNV_OFFSET_32;
    for (std::int32_t v : input) {
        hash = fnv_1a_int32(hash, v);
    }
    return hash;
}

enum class fhash_key_kind { int32_scalar, int32_vector, string };

struct fhash_key_t {
    fhash_key_kind kind = fhash_key_kind::int32_vector;
    std::int32_t int_scalar = 0;
    std::vector<std::int32_t> int_vector;
    std::string str;

    std::uint64_t hash() const {
        switch (kind) {
        case fhash_key_kind::int32_scalar:
            return fnv_1a_int32(FNV_OFFSET_32, int_scalar);
        case fhash_key_kind::int32_vector:
            return fnv_1a_int32_1d(int_vector);
        case fhash_key_kind::string:
            return fnv_1a_string(FNV_OFFSET_32, str);
        }
        return FNV_OFFSET_32;
    }

    bool equals(const fhash_key_t& other) const {
        if (kind != other.kind) {
            return false;
        }
        switch (kind) {
        case fhash_key_kind::int32_scalar:
            return int_scalar == other.int_scalar;
        case fhash_key_kind::int32_vector:
            return int_vector == other.int_vector;
        case fhash_key_kind::string:
            return str == other.str;
        }
        return false;
    }
};

inline fhash_key_t key(const std::vector<int>& source) {
    fhash_key_t k;
    k.kind = fhash_key_kind::int32_vector;
    k.int_vector.reserve(source.size());
    for (int v : source) {
        k.int_vector.push_back(static_cast<std::int32_t>(v));
    }
    return k;
}

inline fhash_key_t key(int source) {
    fhash_key_t k;
    k.kind = fhash_key_kind::int32_scalar;
    k.int_scalar = static_cast<std::int32_t>(source);
    return k;
}

inline fhash_key_t key(const std::string& source) {
    fhash_key_t k;
    k.kind = fhash_key_kind::string;
    k.str = source;
    return k;
}

class fhash_tbl_t {
public:
    fhash_tbl_t() = default;

    void allocate(int size = FHASH_DEFAULT_ALLOCATION) {
        buckets_.assign(static_cast<std::size_t>(size), {});
        allocated_ = true;
    }

    void set(const fhash_key_t& key, const std::any& value) {
        if (!allocated_) {
            allocate();
        }
        const std::size_t index = bucket_index(key);
        auto& chain = buckets_[index];
        for (auto& entry : chain) {
            if (entry.first.equals(key)) {
                entry.second = value;
                return;
            }
        }
        chain.emplace_back(key, value);
    }

    void get_raw(const fhash_key_t& key, std::any& value, int* stat = nullptr) const {
        if (stat) {
            *stat = 0;
        }
        if (!allocated_) {
            if (stat) {
                *stat = FHASH_EMPTY_TABLE;
            }
            return;
        }
        const std::size_t index = bucket_index(key);
        for (const auto& entry : buckets_[index]) {
            if (entry.first.equals(key)) {
                value = entry.second;
                return;
            }
        }
        if (stat) {
            *stat = FHASH_KEY_NOT_FOUND;
        }
    }

    bool hasKey(const fhash_key_t& key) const {
        int stat = 0;
        check_key(key, stat);
        return stat == 0;
    }

    void check_key(const fhash_key_t& key, int& stat) const {
        stat = 0;
        if (!allocated_) {
            stat = FHASH_EMPTY_TABLE;
            return;
        }
        const std::size_t index = bucket_index(key);
        for (const auto& entry : buckets_[index]) {
            if (entry.first.equals(key)) {
                return;
            }
        }
        stat = FHASH_KEY_NOT_FOUND;
    }

    void unset(const fhash_key_t& key, int* stat = nullptr) {
        if (stat) {
            *stat = 0;
        }
        if (!allocated_) {
            if (stat) {
                *stat = FHASH_EMPTY_TABLE;
            }
            return;
        }
        const std::size_t index = bucket_index(key);
        auto& chain = buckets_[index];
        for (auto it = chain.begin(); it != chain.end(); ++it) {
            if (it->first.equals(key)) {
                chain.erase(it);
                return;
            }
        }
        if (stat) {
            *stat = FHASH_KEY_NOT_FOUND;
        }
    }

protected:
    std::vector<std::vector<std::pair<fhash_key_t, std::any>>> buckets_;
    bool allocated_ = false;

    std::size_t bucket_index(const fhash_key_t& key) const {
        const std::uint64_t h = key.hash();
        return static_cast<std::size_t>(h % static_cast<std::uint64_t>(buckets_.size()));
    }
};

} // namespace fhash_m

#endif // FHASH_M_H
