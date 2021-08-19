#ifndef JULY_C2_SYMMETRIZER_H
#define JULY_C2_SYMMETRIZER_H

#include "Task.h"
#include "cmath"

Task& c2_symmetrizer(Task& T, int pairs) {

    if (T.is_C2_symmetrized) {
        return T;
    }

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

    int max_repr = 2;

    auto symmetrized_lex = [&cumulative_product, &T, pairs](const unsigned long lex){
        std::vector<int> nzs(T.mults.size());
        for (int i = 0; i < T.mults.size(); ++i) {
            nzs[i] = (lex % cumulative_product[i]) / cumulative_product[i + 1];
        }
        for (int k = 0; k < pairs; ++k) {
            std::swap(nzs[2 * k], nzs[2 * k + 1]);
        }
        unsigned long symm_lex = 0;
        for (int i = 0; i < T.mults.size(); ++i) {
            symm_lex += nzs[i] * cumulative_product[i + 1];
        }
        return symm_lex;
    };

    auto projector = [&symmetrized_lex, max_repr](std::map<Index, Coefficient> & m){
        std::vector<std::map<Index, Coefficient>> projections(max_repr);
        std::map<Index, Coefficient> m_symm;
        for (auto& p : m) {
            unsigned long symm_lex = symmetrized_lex(p.first);
            m_symm[symm_lex] = p.second;
        }
        auto it_m = m.begin();
        auto it_m_symm = m_symm.begin();
        while (it_m != m.end()) {
            Coefficient old_coeff = it_m->second;
            projections[0][it_m->first] += 1 * old_coeff;
            projections[0][it_m_symm->first] += 1 * old_coeff;
            projections[1][it_m->first] += 1 * old_coeff;
            projections[1][it_m_symm->first] += -1 * old_coeff;
            ++it_m;
            ++it_m_symm;
        }
        for (std::map<Index, Coefficient>& mm : projections) {
            for (auto i = mm.begin(), last = mm.end(); i != last; ) {
                if (std::abs(i->second) < 0.001) {
                    i = mm.erase(i);
                } else {
                    ++i;
                }
            }
        }
        return projections;
    };

    while (T.blocks.front().representation == -1) {
        Subspace& Ss_front = T.blocks.front();
        Subspace Ss_blank = T.blocks.front();
        Ss_blank.basis.clear();
        std::vector<int> repr_to_block(max_repr, -1);

        for (int i = 0; i < Ss_front.basis.size(); ++i) {
            std::vector<std::map<Index, Coefficient>> projections = projector(Ss_front.basis[i]);
            for (int repr = 0; repr < max_repr; ++repr) {
                if (!projections[repr].empty()) {
                    if (repr_to_block[repr] == -1) {
                        T.blocks.push_back(Ss_blank);
                        T.blocks.back().representation = repr;
                        repr_to_block[repr] = T.blocks.size() - 1;
                    }
                    int j = repr_to_block[repr];
                    T.blocks[j].basis.emplace_back(std::move(projections[repr]));
                }
            }
        }

        T.blocks.pop_front();
    }

    T.is_C2_symmetrized = true;
    return T;
}


#endif //JULY_C2_SYMMETRIZER_H
