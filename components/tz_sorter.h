#ifndef JULY_TZ_SORTER_H
#define JULY_TZ_SORTER_H

#include "Task.h"
#include <numeric>
#include <unordered_map>

Task& tz_sorter(Task& T) {

    if (T.is_Tz_sorted) {
        return T;
    }

    int max_ntz_proj = std::accumulate(T.mults.begin(), T.mults.end(), 1 - T.mults.size());

    int tensor_size = 1;
    for (int mult : T.mults) {
        tensor_size *= mult;
    }

    std::vector<int> cumulative_product;
    cumulative_product.resize(T.mults.size() + 1, 0);
    cumulative_product[0] = tensor_size;
    for (int k = 1; k < T.mults.size() + 1; ++k) {
        cumulative_product[k] = cumulative_product[k - 1] / T.mults[k - 1];
    }

    auto lex_to_ntzproj = [&cumulative_product, &T](unsigned long lex) {
        int ntz_proj = 0;
        for (int i = 0; i < T.mults.size(); ++i) {
            ntz_proj += (lex % cumulative_product[i]) / cumulative_product[i + 1];
        }
        return ntz_proj;
    };

    while (T.blocks.front().n_proj == -1) {
        Subspace& Ss_front = T.blocks.front();
        Subspace Ss_blank = T.blocks.front();
        Ss_blank.basis.clear();
        std::vector<int> ntz_proj_to_block(max_ntz_proj, -1);

        for (int i = 0; i < Ss_front.basis.size(); ++i) {
            int ntz_proj = lex_to_ntzproj(Ss_front.basis[i].begin()->first);
//            int ntz_proj = lex_to_ntzproj(Ss_front.basis[i][0].index);
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

#endif //JULY_TZ_SORTER_H
