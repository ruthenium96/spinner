#include "c2_symmetrizer.h"

Symmetrizer::Symmetrizer(const Spaces::Indexes& indexes, int pairs_)
    : indexes_(indexes), pairs(pairs_) {
    // TODO: copy it from Group object:
    max_repr = 2;
}

Space& Symmetrizer::operator()(Space& space) {
    if (space.is_C2_symmetrized) {
        return space;
    }

    while (space.blocks.front().representation == -1) {
        Subspace& subspace_parent = space.blocks.front();
        Subspace subspace_child = space.blocks.front();
        subspace_child.basis.clear();
        std::vector<int> repr_to_block(max_repr, -1);

        std::unordered_map<size_t, size_t> visited;

        for (int i = 0; i < subspace_parent.basis.size(); ++i) {
            if (!is_in_hash_table(subspace_parent.basis[i], visited)) {
                std::vector<Decomposition> projections =
                    projector(subspace_parent.basis[i], visited);
                for (int repr = 0; repr < max_repr; ++repr) {
                    if (!projections[repr].empty()) {
                        if (repr_to_block[repr] == -1) {
                            space.blocks.push_back(subspace_child);
                            space.blocks.back().representation = repr;
                            repr_to_block[repr] = space.blocks.size() - 1;
                        }
                        int j = repr_to_block[repr];
                        space.blocks[j].basis.emplace_back(std::move(projections[repr]));
                    }
                }
            }
        }

        space.blocks.pop_front();
    }

    space.is_C2_symmetrized = true;
    return space;
}

unsigned long Symmetrizer::symmetrized_lex(const unsigned long lex) const {
    std::vector<int> nzs = indexes_.lex_to_nzs(lex);
    for (int k = 0; k < pairs; ++k) {
        std::swap(nzs[2 * k], nzs[2 * k + 1]);
    }
    unsigned long symm_lex = indexes_.nzs_to_lex(nzs);
    return symm_lex;
};

std::vector<Decomposition> Symmetrizer::projector(Decomposition & m,
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

void Symmetrizer::add_to_hash_table(Decomposition & m,
                                    std::unordered_map<size_t, size_t>& hs) {
    size_t key_seed = 0;
    size_t value_seed = 11;
    boost::hash_range(key_seed, m.begin(), m.end());
    boost::hash_range(value_seed, m.begin(), m.end());
    hs[key_seed] = value_seed;
}

void Symmetrizer::erase_if_zero(std::vector<Decomposition>& projections) {
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

bool Symmetrizer::is_in_hash_table(const Decomposition & m,
                                   std::unordered_map<size_t, size_t>& hs) {
    size_t key_seed = 0;
    size_t value_seed = 11;
    boost::hash_range(key_seed, m.begin(), m.end());
    boost::hash_range(value_seed, m.begin(), m.end());
    return hs[key_seed] == value_seed;
}