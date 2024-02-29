#ifndef SPINNER_MUMXHASH_H
#define SPINNER_MUMXHASH_H

#include <cstddef>
#include <cstdint>

namespace ankerl {

// We took this hash and mixer from:
// https://github.com/martinus/map_benchmark/blob/master/src/app/mixer.h
// https://github.com/martinus/map_benchmark/blob/master/src/hashes/mumx_hash/Hash.h
inline uint64_t umul128(uint64_t a, uint64_t b, uint64_t* high) noexcept {
#if defined(__SIZEOF_INT128__)
    #if defined(__GNUC__) || defined(__clang__)
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wpedantic"
    using uint128_t = unsigned __int128;
        #pragma GCC diagnostic pop
    #endif

    auto result = static_cast<uint128_t>(a) * static_cast<uint128_t>(b);
    *high = static_cast<uint64_t>(result >> 64U);
    return static_cast<uint64_t>(result);

#elif (defined(_MSC_VER) && SIZE_MAX == UINT64_MAX)
    #include <intrin.h>  // for __umulh
    #pragma intrinsic(__umulh)
    #ifndef _M_ARM64
        #pragma intrinsic(_umul128)
    #endif
    #ifdef _M_ARM64
    *high = __umulh(a, b);
    return ((uint64_t)(a)) * (b);
    #else
    return _umul128(a, b, high);
    #endif
#else
    #error No hardware umul
#endif
}

// 128bit multiply a and b, xor high and low result
inline uint64_t mumx(uint64_t a, uint64_t b) {
    uint64_t h;
    uint64_t l = umul128(a, b, &h);
    return h ^ l;
}

template<typename T>
struct MumxHash {
    size_t operator()(T const& v) const {
        static constexpr auto a = UINT64_C(0x2ca7aea0ebd71d49);
        return ankerl::mumx(v, a);
    }
};

}  // namespace ankerl

#endif  //SPINNER_MUMXHASH_H
