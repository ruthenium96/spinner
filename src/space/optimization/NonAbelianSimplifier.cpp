#include "NonAbelianSimplifier.h"

#include <cassert>
#include <utility>

namespace space::optimization {

Space NonAbelianSimplifier::apply(Space&& space) const {
    assert(!space.getBlocks().empty());
    auto totalSpaceSize = space.getBlocks()[0].decomposition->size_rows();
    std::vector<Subspace> vector_result;
    vector_result.reserve(space.getBlocks().size());
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        vector_result.emplace_back(factories_.createSparseSemiunitaryMatrix(0, totalSpaceSize));
    }

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace_parent = space.getBlocks()[i];
        uint32_t old_dimensionality = subspace_parent.properties.dimensionality;
        if (old_dimensionality == 1) {
            // it is senseless to modify subspaces with old_dimensionality == 1:
            vector_result[i] = std::move(subspace_parent);
            continue;
        }

        BlockProperties block_properties = subspace_parent.properties;
        block_properties.dimensionality = 1;
        block_properties.degeneracy *= old_dimensionality;
        vector_result[i].properties = block_properties;

        for (size_t j = 0; j < subspace_parent.decomposition->size_cols();
             j = j + old_dimensionality) {
            vector_result[i].decomposition->move_vector_from(j, subspace_parent.decomposition);
        }
        subspace_parent.decomposition->clear();
    }
    return Space(std::move(vector_result));
}
NonAbelianSimplifier::NonAbelianSimplifier(quantum::linear_algebra::FactoriesList factories) :
    factories_(std::move(factories)) {}
}  // namespace space::optimization