#include "TzSorter.h"

#include <utility>

Space TzSorter::apply(Space& space) const {
    // It does not make any sense to use tz_sorter twice.
    if (space.history.isTzSorted) {
        return std::move(space);
    }

    std::vector<Subspace> vector_result;
    vector_result.resize(space.blocks.size() * (max_ntz_proj + 1));

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.blocks.size(); ++i) {
        Subspace& subspace_parent = space.blocks[i];

        for (size_t ntz_proj = 0; ntz_proj < max_ntz_proj + 1; ++ntz_proj) {
            BlockProperties block_properties = subspace_parent.properties;
            block_properties.n_proj = ntz_proj;
            vector_result[(max_ntz_proj + 1) * i + ntz_proj].properties = block_properties;
        }

        for (auto & basi : subspace_parent.basis) {
            // Value of total projection is calculated from the first index of map.
            // NB: there is no validation of the fact, that all indexes of decomposition
            // correspond to the same projection value, user should check it yourself.
            uint8_t ntz_proj = converter_.convert_lex_index_to_tz_projection(basi.begin()->first);
            size_t j = (max_ntz_proj + 1) * i + ntz_proj;
            vector_result[j].basis.emplace_back(std::move(basi));
        }

        subspace_parent.basis.clear();
    }

    return Space(std::move(vector_result));
}

TzSorter::TzSorter(spaces::LexicographicIndexConverter indexes) : converter_(std::move(indexes)) {
    // We want to get 2T + 1 (projections are counted from zero to multiplicity),
    // where T = sum_{1}^{N} S_i. So 2T + 1 = sum_{1}^{N} (2S_i + 1) - N + 1.
    max_ntz_proj = std::accumulate(converter_.mults_.begin(),
                                   converter_.mults_.end(),
                                   1 - converter_.mults_.size());
}
