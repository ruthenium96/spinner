#ifndef JULY_SYMMETRIZER_H
#define JULY_SYMMETRIZER_H

#include <boost/functional/hash.hpp>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "src/group/Group.h"
#include "src/space/Space.h"

namespace space::optimization {

class Symmetrizer {
  public:
    Symmetrizer(lexicographic::IndexConverter converter, group::Group group);
    Space apply(Space&& space) const;

    std::vector<UnitarySparseMatrix>
    get_symmetrical_projected_decompositions(Subspace& subspace, uint32_t index_of_vector) const;

    // TODO: these functions are not about symmetrization.
    //  Should we refactor them and create a new class?
    static void increment_visited(
        const UnitarySparseMatrix& decomposition,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, uint8_t>& hs);

    static uint8_t count_how_many_orbit_was_visited(
        const UnitarySparseMatrix& decomposition,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, uint8_t>& hs);

    static bool is_orthogonal_to_others(
        const UnitarySparseMatrix& decomposition_from,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, std::vector<size_t>>& hs,
        const Subspace& subspace_to);
    static void move_vector_and_remember_it(
        UnitarySparseMatrix& decomposition_from,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, std::vector<size_t>>& hs,
        Subspace& subspace_to);

  private:
    const lexicographic::IndexConverter converter_;
    const group::Group group_;
};
}  // namespace space::optimization

#endif  // JULY_SYMMETRIZER_H
