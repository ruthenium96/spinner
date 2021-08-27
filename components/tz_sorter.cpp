#include "tz_sorter.h"

Space &Tz_Sorter::operator()(Space &space) const {
    // It does not make any sense to use tz_sorter twice.
    if (space.is_Tz_sorted) {
        return space;
    }
    while (space.blocks.front().n_proj == -1) {
        Subspace& subspace_parent = space.blocks.front();
        Subspace subspace_child = space.blocks.front();
        subspace_child.basis.clear();
        // it is a mapping between Subspace with specific projection value and position into deque
        std::vector<int> ntz_proj_to_block(max_ntz_proj, -1);

        for (int i = 0; i < subspace_parent.basis.size(); ++i) {
            // Value of total projection is calculated from the first index of map.
            // NB: there is no validation of the fact, that all indexes of decomposition
            // correspond to the same projection value, user should check it yourself.
            int ntz_proj = indexes_.lex_to_ntzproj(subspace_parent.basis[i].begin()->first);
            // if it is the first basis vector with this ntz_proj, create new Subspace in deque
            if (ntz_proj_to_block[ntz_proj] == -1) {
                space.blocks.push_back(subspace_child);
                space.blocks.back().n_proj = ntz_proj;
                ntz_proj_to_block[ntz_proj] = space.blocks.size() - 1;
            }
            int j = ntz_proj_to_block[ntz_proj];
            space.blocks[j].basis.emplace_back(std::move(subspace_parent.basis[i]));
        }

        space.blocks.pop_front();
    }

    space.is_Tz_sorted = true;
    return space;
}

Tz_Sorter::Tz_Sorter(const Spaces::Indexes& indexes) : indexes_(indexes) {
    // We want to get 2T + 1 (projections are counted from zero to multiplicity), where T = sum_{1}^{N} S_i.
    // So 2T + 1 = sum_{1}^{N} (2S_i + 1) - N + 1
    max_ntz_proj = std::accumulate(indexes.mults_.begin(),
                                   indexes.mults_.end(),
                                   1 - indexes.mults_.size());
}
