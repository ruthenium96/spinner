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

    std::vector<Subspace> get_symmetrical_projected_decompositions(Subspace& subspace, uint32_t index_of_vector) const;

    // TODO: these functions are not about symmetrization.
    //  Should we refactor them and create a new class?
    static void increment_visited(const Subspace& subspace, uint32_t index_of_vector, std::unordered_map<uint32_t , uint8_t>& hs);

    static uint8_t count_how_many_orbit_was_visited(const Subspace& subspace, uint32_t index_of_vector, std::unordered_map<uint32_t , uint8_t>& hs);

    static void add_vector_if_orthogonal_to_others(Subspace& subspace_from, uint32_t index_of_vector, std::unordered_map<uint32_t, std::vector<size_t>>& hs,
                                                    Subspace& subspace_to);


  private:
    const spaces::LexicographicIndexConverter converter_;
    const Group group_;
};

#endif // JULY_SYMMETRIZER_H
