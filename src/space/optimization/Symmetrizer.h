#ifndef SPINNER_SYMMETRIZER_H
#define SPINNER_SYMMETRIZER_H

#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "src/common/index_converter/AbstractIndexConverter.h"
#include "src/entities/data_structures/FactoriesList.h"
#include "src/group/Group.h"
#include "src/space/Space.h"

namespace space::optimization {

class Symmetrizer {
  public:
    Symmetrizer(
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        group::Group group,
        quantum::linear_algebra::FactoriesList factories);
    Space apply(Space&& space) const;

  private:
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>>
    get_symmetrical_projected_decompositions(Subspace& subspace, uint32_t index_of_vector) const;

    // TODO: these functions are not about symmetrization.
    //  Should we refactor them and create a new class?
    static void increment_visited(
        const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
            decomposition,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, uint8_t>& hs);

    static uint8_t count_how_many_orbit_was_visited(
        const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
            decomposition,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, uint8_t>& hs);

    static bool is_orthogonal_to_others(
        const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
            decomposition_from,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, std::vector<size_t>>& hs,
        const Subspace& subspace_to);
    static void move_vector_and_remember_it(
        std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
            decomposition_from,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, std::vector<size_t>>& hs,
        Subspace& subspace_to);

    std::shared_ptr<const index_converter::AbstractIndexConverter> converter_;
    const group::Group group_;
    const quantum::linear_algebra::FactoriesList factories_;
};
}  // namespace space::optimization

#endif  // SPINNER_SYMMETRIZER_H
