#include "ClebshGordanCalculator.h"

#include <omp.h>
#include <wignerSymbols.h>

namespace spin_algebra {

ClebshGordanCalculator::ClebshGordanCalculator() {
    hashed_CGs_for_all_threads = std::vector<Map>(omp_get_max_threads());
}

double ClebshGordanCalculator::clebsh_gordan_coefficient(
    double l1,
    double l2,
    double l3,
    double m1,
    double m2) const {
    if (l1 + l2 < l3 || std::abs(l1 - l2) > l3) {
        return 0;
    }
    auto key = to_key(l1, l2, l3, m1, m2);
    // todo: calculate omp_get_thread_num() outside of this function?
    Map& hashed_CGs_per_thread = hashed_CGs_for_all_threads.at(omp_get_thread_num());

    if (hashed_CGs_per_thread.contains(key)) {
        return hashed_CGs_per_thread.at(key);
    } else {
        double value = WignerSymbols::clebschGordan(l1, l2, l3, m1, m2, m1 + m2);
        hashed_CGs_per_thread.emplace(key, value);
        return value;
    }
}

uint64_t
ClebshGordanCalculator::to_key(double l1, double l2, double l3, double m1, double m2) noexcept {
    // NB: each value cannot be bigger than BASE
    const uint64_t BASE = 4096;
    auto ul1 = (uint64_t)(2 * l1);
    auto ul2 = (uint64_t)(2 * l2);
    auto ul3 = (uint64_t)(2 * l3);
    auto um1 = (uint64_t)(m1 + l1);
    auto um2 = (uint64_t)(m2 + l2);
    uint64_t answer = 0;
    answer |= ul1;
    answer |= ul2 * BASE;
    answer |= ul3 * BASE * BASE;
    answer |= um1 * BASE * BASE * BASE;
    answer |= um2 * BASE * BASE * BASE * BASE;
    return answer;
}
}  // namespace spin_algebra