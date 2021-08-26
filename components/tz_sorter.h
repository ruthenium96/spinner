#ifndef JULY_TZ_SORTER_H
#define JULY_TZ_SORTER_H

#include "common/Space.h"
#include <numeric>
#include <utility>

class Tz_Sorter {
public:
    explicit Tz_Sorter(std::vector<int> mults_) : mults(std::move(mults_)) {
        // We want to get 2T + 1 (projections are counted from zero to multiplicity), where T = sum_{1}^{N} S_i.
        // So 2T + 1 = sum_{1}^{N} (2S_i + 1) - N + 1
        max_ntz_proj = std::accumulate(mults.begin(), mults.end(), 1 - mults.size());

        // It is the size of our projection-based space.
        tensor_size = 1;
        for (int mult : mults) {
            tensor_size *= mult;
        }

        cumulative_product.resize(mults.size() + 1, 0);
        cumulative_product[0] = tensor_size;
        for (int k = 1; k < mults.size() + 1; ++k) {
            cumulative_product[k] = cumulative_product[k - 1] / mults[k - 1];
        }
    }

    int lex_to_ntzproj (unsigned long lex) const {
        int ntz_proj = 0;
        for (int i = 0; i < mults.size(); ++i) {
            ntz_proj += (lex % cumulative_product[i]) / cumulative_product[i + 1];
        }
        return ntz_proj;
    };

    Space& operator()(Space& space) const {
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
                int ntz_proj = lex_to_ntzproj(subspace_parent.basis[i].begin()->first);
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

private:
    const std::vector<int> mults;
    std::vector<int> cumulative_product;
    int max_ntz_proj;
    int tensor_size;
};

#endif //JULY_TZ_SORTER_H
