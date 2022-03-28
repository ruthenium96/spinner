#include "TzSorter.h"

#include <utility>

namespace space::optimization {

Space TzSorter::apply(Space&& space) const {
    std::vector<Subspace> vector_result;
    vector_result.resize(space.blocks.size() * (max_ntz_proj + 1));

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.blocks.size(); ++i) {
        Subspace& subspace_parent = space.blocks[i];

        for (size_t ntz_proj = 0; ntz_proj < max_ntz_proj + 1; ++ntz_proj) {
            BlockProperties block_properties = subspace_parent.properties;
            block_properties.n_proj = ntz_proj;
            vector_result[(max_ntz_proj + 1) * i + ntz_proj].properties = block_properties;
            vector_result[(max_ntz_proj + 1) * i + ntz_proj].decomposition.tensor_size =
                subspace_parent.decomposition.tensor_size;
        }

        for (uint32_t l = 0; l < subspace_parent.decomposition.size(); ++l) {
            // Value of total projection is calculated from the first index of map.
            // NB: there is no validation of the fact, that all indexes of decomposition
            // correspond to the same projection value, user should check it yourself.
            uint32_t index = subspace_parent.decomposition.GetNewIterator(l)->getNext().index;
            uint8_t ntz_proj = converter_.convert_lex_index_to_tz_projection(index);
            size_t j = (max_ntz_proj + 1) * i + ntz_proj;
            vector_result[j].decomposition.move_vector_from(l, subspace_parent.decomposition);
        }

        subspace_parent.decomposition.clear();
    }

    return Space(std::move(vector_result));
}

TzSorter::TzSorter(lexicographic::IndexConverter indexes) : converter_(std::move(indexes)) {
    max_ntz_proj = converter_.get_max_ntz_proj();
}
}  // namespace space::optimization