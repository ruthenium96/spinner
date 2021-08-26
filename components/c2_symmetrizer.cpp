
#include "c2_symmetrizer.h"

Symmetrizer::Symmetrizer(std::vector<int> mults_, int pairs_)
    : mults(mults_), pairs(pairs_) {
    tensor_size = 1;
    for (int mult : mults) {
        tensor_size *= mult;
    }

    cumulative_product.resize(mults.size() + 1, 0);
    cumulative_product[0] = tensor_size;
    for (int k = 1; k < mults.size() + 1; ++k) {
        cumulative_product[k] = cumulative_product[k - 1] / mults[k - 1];
    }

    // TODO: copy it from Group object:
    max_repr = 2;
}

Task& Symmetrizer::operator()(Task& T) {
    if (T.is_C2_symmetrized) {
        return T;
    }

    while (T.blocks.front().representation == -1) {
        Subspace& Ss_front = T.blocks.front();
        Subspace Ss_blank = T.blocks.front();
        Ss_blank.basis.clear();
        std::vector<int> repr_to_block(max_repr, -1);

        std::unordered_map<size_t, size_t> visited;

        for (int i = 0; i < Ss_front.basis.size(); ++i) {
            if (!is_in_hash_table(Ss_front.basis[i], visited)) {
                std::vector<std::map<Index, Coefficient>> projections =
                    projector(Ss_front.basis[i], visited);
                for (int repr = 0; repr < max_repr; ++repr) {
                    if (!projections[repr].empty()) {
                        if (repr_to_block[repr] == -1) {
                            T.blocks.push_back(Ss_blank);
                            T.blocks.back().representation = repr;
                            repr_to_block[repr] = T.blocks.size() - 1;
                        }
                        int j = repr_to_block[repr];
                        T.blocks[j].basis.emplace_back(
                            std::move(projections[repr]));
                    }
                }
            }
        }

        T.blocks.pop_front();
    }

    T.is_C2_symmetrized = true;
    return T;
}

unsigned long Symmetrizer::symmetrized_lex(const unsigned long lex) const {
    std::vector<int> nzs(mults.size());
    for (int i = 0; i < mults.size(); ++i) {
        nzs[i] = (lex % cumulative_product[i]) / cumulative_product[i + 1];
    }
    for (int k = 0; k < pairs; ++k) {
        std::swap(nzs[2 * k], nzs[2 * k + 1]);
    }
    unsigned long symm_lex = 0;
    for (int i = 0; i < mults.size(); ++i) {
        symm_lex += nzs[i] * cumulative_product[i + 1];
    }
    return symm_lex;
};

std::vector<std::map<Index, Coefficient>>
Symmetrizer::projector(std::map<Index, Coefficient>& m,
                       std::unordered_map<size_t, size_t>& hs) {
    std::vector<std::map<Index, Coefficient>> projections(max_repr);
    std::map<Index, Coefficient> m_symm;
    for (auto& p : m) {
        unsigned long symm_lex = symmetrized_lex(p.first);
        m_symm[symm_lex] = p.second;
    }

    add_to_hash_table(m, hs);
    add_to_hash_table(m_symm, hs);

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
    erase_if_zero(projections);
    return projections;
};

void
Symmetrizer::add_to_hash_table(std::map<Index, Coefficient>& m,
                               std::unordered_map<size_t, size_t>& hs) {
    size_t key_seed = 0;
    size_t value_seed = 11;
    boost::hash_range(key_seed, m.begin(), m.end());
    boost::hash_range(value_seed, m.begin(), m.end());
    hs[key_seed] = value_seed;
}

void Symmetrizer::erase_if_zero(
    std::vector<std::map<Index, Coefficient>>& projections) {
    for (std::map<Index, Coefficient>& mm : projections) {
        for (auto i = mm.begin(), last = mm.end(); i != last;) {
            if (std::abs(i->second) < 0.001) {
                i = mm.erase(i);
            } else {
                ++i;
            }
        }
    }
}

bool
Symmetrizer::is_in_hash_table(const std::map<Index, Coefficient>& m,
                              std::unordered_map<size_t, size_t>& hs) {
    size_t key_seed = 0;
    size_t value_seed = 11;
    boost::hash_range(key_seed, m.begin(), m.end());
    boost::hash_range(value_seed, m.begin(), m.end());
    return hs[key_seed] == value_seed;
}
