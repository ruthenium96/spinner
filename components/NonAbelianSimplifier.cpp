#include "NonAbelianSimplifier.h"

Space NonAbelianSimplifier::apply(Space &space) const {
    std::vector<Subspace> vector_result;
    vector_result.resize(space.blocks.size());
    Space::History history_result = space.history;

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.blocks.size(); ++i) {
        Subspace& subspace_parent = space.blocks[i];
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

        for (size_t j = 0; j < subspace_parent.basis.size(); j = j + old_dimensionality) {
            vector_result[i].basis.emplace_back(std::move(subspace_parent.basis[j]));
        }
        subspace_parent.basis.clear();
    }
    return Space(std::move(vector_result), history_result);
}