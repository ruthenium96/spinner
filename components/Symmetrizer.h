#ifndef JULY_SYMMETRIZER_H
#define JULY_SYMMETRIZER_H

#include "cmath"
#include "entities/Space.h"
#include "unordered_map"
#include "unordered_set"
#include <boost/functional/hash.hpp>
#include <vector>
#include <groups/Group.h>

class Symmetrizer {
  public:
    Symmetrizer(spaces::LexicographicIndexConverter converter, Group group);
    Space apply(Space& space) const;

    std::vector<std::vector<DecompositionMap>> get_symmetrical_projected_decompositions(DecompositionMap& m) const;

    // TODO: these functions are not about symmetrization.
    //  Should we refactor them and create a new class?
    static void increment_in_hash_table(DecompositionMap& m, std::unordered_map<uint32_t , uint8_t>& hs);

    static void erase_if_zero(std::vector<std::vector<DecompositionMap>>& projections);

    static uint8_t count_in_hash_table(const DecompositionMap& m, std::unordered_map<uint32_t , uint8_t>& hs);

    static void add_vector_if_orthogonal_to_others(DecompositionMap& m, std::unordered_map<uint32_t, std::vector<size_t>>& hs,
                                                    std::vector<DecompositionMap>& basis);


  private:
    const spaces::LexicographicIndexConverter converter_;
    const Group group_;
};

#endif // JULY_SYMMETRIZER_H
