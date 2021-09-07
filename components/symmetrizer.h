#ifndef JULY_SYMMETRIZER_H
#define JULY_SYMMETRIZER_H

#include "cmath"
#include "entities/Space.h"
#include "unordered_map"
#include <boost/functional/hash.hpp>
#include <vector>
#include <groups/Group.h>

class Symmetrizer {
  public:
    Symmetrizer(const spaces::LexicographicIndexWorker& indexes, const Group& group);
    Space& operator()(Space& space) const;

    std::vector<std::vector<Decomposition>> get_symmetrical_projected_decompositions(Decomposition& m) const;

    // TODO: these functions are not about symmetrization.
    //  Should we refactor them and create a new class?
    static void increment_in_hash_table(Decomposition& m, std::unordered_map<uint32_t , uint8_t>& hs);

    static void erase_if_zero(std::vector<std::vector<Decomposition>>& projections);

    static uint8_t count_in_hash_table(const Decomposition& m, std::unordered_map<uint32_t , uint8_t>& hs);

  private:
    const spaces::LexicographicIndexWorker& indexes_;
    const Group& group_;
};

#endif // JULY_SYMMETRIZER_H
