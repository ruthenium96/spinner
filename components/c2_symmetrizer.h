#ifndef JULY_C2_SYMMETRIZER_H
#define JULY_C2_SYMMETRIZER_H

#include "Task.h"
#include "cmath"
#include "unordered_set"
#include "unordered_map"
#include <boost/functional/hash.hpp>

struct MapHash {
    std::size_t operator()(const std::map<Index, Coefficient>& m) const noexcept
    {
        return boost::hash_range(m.begin(), m.end());
    }
};

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

    auto projector = [&symmetrized_lex, max_repr](std::map<Index, Coefficient> & m,
            std::unordered_map<size_t, size_t> & hs){
//            std::unordered_set<std::map<Index, Coefficient>, MapHash> & hs){
        std::vector<std::map<Index, Coefficient>> projections(max_repr);
        std::map<Index, Coefficient> m_symm;
        for (auto& p : m) {
            unsigned long symm_lex = symmetrized_lex(p.first);
            m_symm[symm_lex] = p.second;
        }
        size_t key_seed = 0;
        size_t value_seed = 11;
        boost::hash_range(key_seed, m.begin(), m.end());
        boost::hash_range(value_seed, m.begin(), m.end());
        hs[key_seed] = value_seed;

        key_seed = 0;
        value_seed = 11;
        boost::hash_range(key_seed, m_symm.begin(), m_symm.end());
        boost::hash_range(value_seed, m_symm.begin(), m_symm.end());
        hs[key_seed] = value_seed;

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

//        std::unordered_set<std::map<Index, Coefficient>, MapHash> visited;
        std::unordered_map<size_t, size_t> visited;

        for (int i = 0; i < Ss_front.basis.size(); ++i) {
            size_t key_seed = 0;
            size_t value_seed = 11;
            boost::hash_range(key_seed, Ss_front.basis[i].begin(), Ss_front.basis[i].end());
            boost::hash_range(value_seed, Ss_front.basis[i].begin(), Ss_front.basis[i].end());

            if (visited[key_seed] != value_seed) {
                std::vector<std::map<Index, Coefficient>> projections = projector(Ss_front.basis[i], visited);
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
        }

        T.blocks.pop_front();
    }

    T.is_C2_symmetrized = true;
    return T;
}


#endif //JULY_C2_SYMMETRIZER_H
