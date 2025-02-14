#include "NonMinimalProjectionsEliminator.h"

namespace space::optimization {

NonMinimalProjectionsEliminator::NonMinimalProjectionsEliminator(uint32_t max_total_mult) :
	max_total_mult_(max_total_mult) {}

Space NonMinimalProjectionsEliminator::apply(Space&& space) const {
    std::vector<Subspace> vector_result;

    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace_parent = space.getBlocks()[i];
        uint32_t ntz_projection = subspace_parent.properties.n_proj.value();
		auto total_mult = subspace_parent.properties.total_mult.value();
        if ((max_total_mult_ - total_mult) / 2 != ntz_projection) {
            continue;
        }

        subspace_parent.properties.degeneracy *= total_mult;
    
        vector_result.emplace_back(std::move(subspace_parent));
    }
    return Space(std::move(vector_result));

}

}