#include "PositiveProjectionsEliminator.h"

namespace space::optimization {

PositiveProjectionsEliminator::PositiveProjectionsEliminator(uint32_t max_ntz_proj) :
    max_ntz_proj_(max_ntz_proj) {}

Space PositiveProjectionsEliminator::apply(Space&& space) const {
    // 0 1 2 3 => 0 1
    // 0 1 2   => 0 1
    // TODO: fix this minus one
    uint32_t half_of_max_projection = (max_ntz_proj_ - 1) / 2;

    std::vector<Subspace> vector_result;

    for (size_t i = 0; i < space.blocks.size(); ++i) {
        Subspace& subspace_parent = space.blocks[i];
        uint32_t ntz_projection = subspace_parent.properties.n_proj.value();
        if (ntz_projection > half_of_max_projection) {
            continue;
        }

        // TODO: also fix this == 0 instead of == 1
        if (ntz_projection != half_of_max_projection || max_ntz_proj_ % 2 == 0) {
            subspace_parent.properties.degeneracy *= 2;
        }
        vector_result.emplace_back(std::move(subspace_parent));
    }
    return Space(std::move(vector_result));
}
}  // namespace space::optimization