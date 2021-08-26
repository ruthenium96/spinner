#ifndef JULY_TZ_SORTER_H
#define JULY_TZ_SORTER_H

#include "common/Task.h"
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

    Task& operator()(Task& T) const {
        // It does not make any sense to use tz_sorter twice.
        if (T.is_Tz_sorted) {
            return T;
        }
        while (T.blocks.front().n_proj == -1) {
            Subspace& Ss_front = T.blocks.front();
            Subspace Ss_blank = T.blocks.front();
            Ss_blank.basis.clear();
            // it is a mapping between Subspace with specific projection value and position into deque
            std::vector<int> ntz_proj_to_block(max_ntz_proj, -1);

            for (int i = 0; i < Ss_front.basis.size(); ++i) {
                // Value of total projection is calculated from the first index of map.
                // NB: there is no validation of the fact, that all indexes of decomposition
                // correspond to the same projection value, user should check it yourself.
                int ntz_proj = lex_to_ntzproj(Ss_front.basis[i].begin()->first);
                // if it is the first basis vector with this ntz_proj, create new Subspace in deque
                if (ntz_proj_to_block[ntz_proj] == -1) {
                    T.blocks.push_back(Ss_blank);
                    T.blocks.back().n_proj = ntz_proj;
                    ntz_proj_to_block[ntz_proj] = T.blocks.size() - 1;
                }
                int j = ntz_proj_to_block[ntz_proj];
                T.blocks[j].basis.emplace_back(std::move(Ss_front.basis[i]));
            }

            T.blocks.pop_front();
        }

        T.is_Tz_sorted = true;
        return T;
    }

private:
    const std::vector<int> mults;
    std::vector<int> cumulative_product;
    int max_ntz_proj;
    int tensor_size;
};

#endif //JULY_TZ_SORTER_H
