#include "NonAbelianSimplifier.h"

Space NonAbelianSimplifier::apply(Space &space) const {
    std::vector<Subspace> vector_result;
    entities::Entity::History history_result = space.history;

    for (Subspace& subspace_parent : space.blocks) {
        int old_dimensionality = subspace_parent.properties.dimensionality;
        if (old_dimensionality == 1) {
            // it is senseless to modify subspaces with old_dimensionality == 1:
            vector_result.emplace_back(subspace_parent);
            continue;
        }
        vector_result.emplace_back();
        vector_result.back().properties = subspace_parent.properties;
        vector_result.back().properties.dimensionality = 1;
        vector_result.back().properties.degeneracy *= old_dimensionality;

        for (size_t i = 0; i < subspace_parent.basis.size(); i = i + old_dimensionality) {
            vector_result.back().basis.emplace_back(std::move(subspace_parent.basis[i]));
        }
        subspace_parent.basis.clear();
    }
    return Space(vector_result, history_result);
}