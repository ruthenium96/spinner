#include "tz_sorter.h"

Space Tz_Sorter::apply(Space& space) const {
    // It does not make any sense to use tz_sorter twice.
    if (space.history.isTzSorted) {
        return std::move(space);
    }

    std::vector<Subspace> vector_result;
    entities::Entity::History history_result = space.history;

    for (Subspace& subspace_parent : space.blocks) {
        // it is a mapping between Subspace with specific projection value and position into deque
        std::vector<size_t> ntz_proj_to_block(max_ntz_proj, -1);

        for (auto & basi : subspace_parent.basis) {
            // Value of total projection is calculated from the first index of map.
            // NB: there is no validation of the fact, that all indexes of decomposition
            // correspond to the same projection value, user should check it yourself.
            uint8_t ntz_proj = indexes_.convert_lex_index_to_tz_projection(basi.begin()->first);
            // if it is the first basis vector with this ntz_proj, create new Subspace in deque
            if (ntz_proj_to_block[ntz_proj] == -1) {
                vector_result.emplace_back();
                vector_result.back().properties = subspace_parent.properties;
                vector_result.back().properties.n_proj = ntz_proj;
                ntz_proj_to_block[ntz_proj] = vector_result.size() - 1;
            }
            size_t j = ntz_proj_to_block[ntz_proj];

            vector_result[j].basis.emplace_back(std::move(basi));
        }

        subspace_parent.basis.clear();
    }

    history_result.isTzSorted = true;
    return Space(vector_result, history_result);
}

Tz_Sorter::Tz_Sorter(const spaces::LexicographicIndexConverter& indexes) : indexes_(indexes) {
    // We want to get 2T + 1 (projections are counted from zero to multiplicity),
    // where T = sum_{1}^{N} S_i. So 2T + 1 = sum_{1}^{N} (2S_i + 1) - N + 1.
    max_ntz_proj = std::accumulate(indexes.mults_.begin(),
                                   indexes.mults_.end(),
                                   1 - indexes.mults_.size());
}
