#ifndef SPINNER_SYMMETRIZER_H
#define SPINNER_SYMMETRIZER_H

#include <cmath>
#include <unordered_map>
#include <vector>

#include "src/common/index_converter/AbstractIndexPermutator.h"
#include "src/entities/data_structures/FactoriesList.h"
#include "src/group/Group.h"
#include "src/space/Space.h"

namespace space::optimization {

class Symmetrizer {
  public:
    Symmetrizer(
        std::shared_ptr<const index_converter::AbstractIndexPermutator> permutator,
        group::Group group,
        quantum::linear_algebra::FactoriesList factories);
    Space apply(Space&& space) const;

  private:
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>>
    get_symmetrical_projected_decompositions(Subspace& subspace, uint32_t index_of_vector) const;

    // TODO: these functions are not about symmetrization.
    //  Should we refactor them and create a new class?
    static void gram_schmidt_orthogonalize(
        std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& decomposition_from,
        const std::unordered_map<uint32_t, std::vector<size_t>>& indexes_to_vectors_map,
        const Subspace& subspace_to);
    static void gram_schmidt_selforthogonalize(
        std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
            decomposition_from);
    static void move_vector_and_remember_it(
        std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
            decomposition_from,
        uint32_t index_of_vector,
        std::unordered_map<uint32_t, std::vector<size_t>>& hs,
        Subspace& subspace_to);

    std::shared_ptr<const index_converter::AbstractIndexPermutator> permutator_;
    const group::Group group_;
    const quantum::linear_algebra::FactoriesList factories_;
};
}  // namespace space::optimization

#endif  // SPINNER_SYMMETRIZER_H
