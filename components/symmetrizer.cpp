#include "symmetrizer.h"

Symmetrizer::Symmetrizer(const spaces::LexicographicIndexWorker& indexes, const Group& group)
    : indexes_(indexes), group_(group) {
}

Space& Symmetrizer::operator()(Space& space) const {
//    if (space.is_C2_symmetrized) {
//        return space;
//    }

    // save old numbers of symmetrizer usage.
    size_t old_number_of_simmetrization = space.blocks.front().representation.size();

    // while we are working with old subspaces:
    while (space.blocks.front().representation.size() == old_number_of_simmetrization) {
        Subspace& subspace_parent = space.blocks.front();
        Subspace subspace_child = space.blocks.front();
        subspace_child.basis.clear();
        std::vector<size_t> repr_to_block(group_.groupInfo.number_of_representations, -1);

        // it is an auxiliary hash table. It helps excludes exact copies of the vectors.
        std::unordered_map<size_t, size_t> visited;

        for (auto & basi : subspace_parent.basis) {
            // when we work with basi, we actually do all work for the orbit of basi,
            // so we add basi and its orbits to visited,
            // because there is no reason to work with them over and over
            if (!is_in_hash_table(basi, visited)) {
                std::vector<Decomposition> projected_basi = get_symmetrical_projected_decompositions(
                        basi, visited);
                for (int repr = 0; repr < group_.groupInfo.number_of_representations; ++repr) {
                    if (!projected_basi[repr].empty()) {
                        if (repr_to_block[repr] == -1) {
                            space.blocks.push_back(subspace_child);
                            space.blocks.back().representation.emplace_back(repr);
                            repr_to_block[repr] = space.blocks.size() - 1;
                        }
                        size_t j = repr_to_block[repr];
                        space.blocks[j].basis.emplace_back(std::move(projected_basi[repr]));
                    }
                }
            }
        }

        space.blocks.pop_front();
    }

//    space.is_C2_symmetrized = true;
    return space;
}

std::vector<Decomposition> Symmetrizer::get_symmetrical_projected_decompositions(Decomposition & m,
                                                                                 std::unordered_map<size_t, size_t>& hs) const {
    std::vector<Decomposition> projections(group_.groupInfo.number_of_representations);
    std::vector<Decomposition> ms_parent(group_.groupInfo.group_size);

    for (auto& p : m) {
        std::vector<uint8_t> nzs = indexes_.lex_to_nzs(p.first);
        std::vector<std::vector<uint8_t>> permutated_vectors = group_.permutate(nzs);

        for (uint32_t i = 0; i < permutated_vectors.size(); ++i) {
            uint32_t permutated_lex = indexes_.nzs_to_lex(permutated_vectors[i]);
            ms_parent[i][permutated_lex] = p.second;
        }
    }

    for (auto& mm : ms_parent) {
        add_to_hash_table(mm, hs);
    }

    std::vector<Decomposition::iterator> ms_parent_iterator(group_.groupInfo.group_size);
    for (uint32_t i = 0; i < group_.groupInfo.group_size; ++i) {
        ms_parent_iterator[i] = ms_parent[i].begin();
    }
    while (ms_parent_iterator[0] != ms_parent[0].end()) {
        double old_coeff = ms_parent_iterator[0]->second;
        for (uint32_t i = 0; i < group_.groupInfo.number_of_representations; ++i) {
            for (uint32_t j = 0; j < group_.groupInfo.group_size; ++j) {
                projections[i][ms_parent_iterator[j]->first] += group_.groupInfo.coefficients_of_projectors[i][j] * old_coeff;
            }
        }
        for (uint32_t i = 0; i < group_.groupInfo.group_size; ++i) {
            ++ms_parent_iterator[i];
        }
    }

    // this function deletes all pairs (key, value) with value = 0.0 from projections
    erase_if_zero(projections);
    return projections;
}

void Symmetrizer::add_to_hash_table(Decomposition & m,
                                    std::unordered_map<size_t, size_t>& hs) {
    // TODO: It is not the fastest way, can we do something else?
    // TODO: Probably, implementing of something like GCD will be also good (but still slow).
    if (m.begin()->second < 0) {
        for (auto& p : m) {
            p.second *= -1;
        }
    }
    size_t key_seed = 0;
    size_t value_seed = 11;
    boost::hash_range(key_seed, m.begin(), m.end());
    boost::hash_range(value_seed, m.begin(), m.end());
    hs[key_seed] = value_seed;
}

void Symmetrizer::erase_if_zero(std::vector<Decomposition>& projections) {
    for (Decomposition& mm : projections) {
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
