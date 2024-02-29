#ifndef SPINNER_CLEBSHGORDANCALCULATOR_H
#define SPINNER_CLEBSHGORDANCALCULATOR_H

#include <hash_table8.hpp>
#include <vector>

#include "src/common/MumxHash.h"

namespace spin_algebra {
class ClebshGordanCalculator {
  public:
    using Map = emhash8::HashMap<uint64_t, double, ankerl::MumxHash<uint64_t>>;

    ClebshGordanCalculator();

    double clebsh_gordan_coefficient(double l1, double l2, double l3, double m1, double m2) const;
    static uint64_t to_key(double l1, double l2, double l3, double m1, double m2) noexcept;

  private:
    // we use caching of CG values to speed up calculations.
    // real concurrent hashmaps are too slow for our task,
    // thus we will use ordinary ones -- one per each thread,
    // this solution already significantly improves performance
    mutable std::vector<Map> hashed_CGs_for_all_threads;
};
}  // namespace spin_algebra

#endif  //SPINNER_CLEBSHGORDANCALCULATOR_H
