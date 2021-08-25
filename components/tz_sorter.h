#ifndef JULY_TZ_SORTER_H
#define JULY_TZ_SORTER_H

#include "Task.h"
#include <numeric>
#include <utility>

class Tz_Sorter {
public:
    explicit Tz_Sorter(std::vector<int> mults_) : mults(std::move(mults_)) {
        max_ntz_proj = std::accumulate(mults.begin(), mults.end(), 1 - mults.size());

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
        if (T.is_Tz_sorted) {
            return T;
        }
        while (T.blocks.front().n_proj == -1) {
            Subspace& Ss_front = T.blocks.front();
            Subspace Ss_blank = T.blocks.front();
            Ss_blank.basis.clear();
            std::vector<int> ntz_proj_to_block(max_ntz_proj, -1);

            for (int i = 0; i < Ss_front.basis.size(); ++i) {
                int ntz_proj = lex_to_ntzproj(Ss_front.basis[i].begin()->first);
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
