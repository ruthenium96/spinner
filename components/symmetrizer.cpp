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
        std::vector<size_t> repr_to_block(group_.info.number_of_representations, -1);

        // it is an auxiliary hash table. It helps to calculate each orbit only once (see below).
        std::unordered_map<uint32_t, uint8_t> visited;

        for (auto & basi : subspace_parent.basis) {
            uint8_t dimension_of_parent;
            if (subspace_parent.representation.empty()) {
                dimension_of_parent = 1;
            } else {
                // TODO: now we can use only the same groups, because we don't have Properties class
                //  here we have to use _the previous_ group, not the current one
                //  another idea: keep dimension (and degeneracy) of subspace
                uint8_t last_representation = subspace_parent.representation.back();
                dimension_of_parent = group_.info.dimension_of_representation[last_representation];
            }
            // when we work with basi, we actually do all work for the orbit of basi,
            // so we add basi and its orbits to visited,
            // because there is no reason to work with them over and over
            if (count_in_hash_table(basi, visited) < dimension_of_parent) {
                std::vector<std::vector<Decomposition>> projected_basi = get_symmetrical_projected_decompositions(basi);
                // TODO: we actually have to add only first (full symmetrical) representation,
                //  but this implementation knows where it is
                increment_in_hash_table(projected_basi[0][0], visited);
                for (uint8_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
                    std::unordered_map<uint32_t, uint8_t> added;
                    uint8_t dimension_of_child = group_.info.dimension_of_representation[repr];
                    for (uint8_t k = 0; k < group_.info.number_of_projectors_of_representation[repr]; ++k) {
                        // check if the Decomposition is empty:
                        if (!projected_basi[repr][k].empty()) {
                            if (repr_to_block[repr] == -1) {
                                space.blocks.push_back(subspace_child);
                                space.blocks.back().representation.emplace_back(repr);
                                repr_to_block[repr] = space.blocks.size() - 1;
                            }
                            size_t j = repr_to_block[repr];
                            if (count_in_hash_table(projected_basi[repr][k], added) < dimension_of_child) {
                                increment_in_hash_table(projected_basi[repr][k], added);
                                space.blocks[j].basis.emplace_back(std::move(projected_basi[repr][k]));
                            }
                        }
                    }
                }
            }
        }

        space.blocks.pop_front();
    }

//    space.is_C2_symmetrized = true;
    return space;
}

std::vector<std::vector<Decomposition>> Symmetrizer::get_symmetrical_projected_decompositions(Decomposition & m) const {
    // it is a set (partitioned by representations) of all projected decompositions:
    std::vector<std::vector<Decomposition>> projections(group_.info.number_of_representations);
    for (uint8_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
        projections[repr].resize(group_.info.number_of_projectors_of_representation[repr]);
    }

    for (auto& p : m) {
        std::vector<uint8_t> nzs = indexes_.lex_to_nzs(p.first);
        std::vector<std::vector<uint8_t>> permutated_vectors = group_.permutate(nzs);

        for (uint8_t g = 0; g < group_.info.group_size; ++g) {
            uint32_t permutated_lex = indexes_.nzs_to_lex(permutated_vectors[g]);
            for (uint8_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
                for (uint8_t projector = 0; projector < group_.info.number_of_projectors_of_representation[repr]; ++projector) {
                    projections[repr][projector][permutated_lex] += group_.info.coefficients_of_projectors[repr][projector][g] * p.second;
                }
            }
        }
    }

    // this function deletes all pairs (key, value) with value = 0.0 from projections
    erase_if_zero(projections);
    return projections;
}

void Symmetrizer::increment_in_hash_table(Decomposition & m,
                                    std::unordered_map<uint32_t , uint8_t>& hs) {
    for (const auto p : m) {
        if (hs.find(p.first) == hs.end()) {
            hs[p.first] = 1;
        } else {
            ++hs[p.first];
        }
    }
}

void Symmetrizer::erase_if_zero(std::vector<std::vector<Decomposition>>& projections) {
    for (auto & v : projections) {
        for (Decomposition& mm : v) {
            for (auto i = mm.begin(), last = mm.end(); i != last;) {
                if (std::abs(i->second) < 0.001) {
                    i = mm.erase(i);
                } else {
                    ++i;
                }
            }
        }
    }
}

uint8_t Symmetrizer::count_in_hash_table(const Decomposition & m,
                                   std::unordered_map<uint32_t , uint8_t>& hs) {
    uint8_t counter = 0;
    for (const auto p : m) {
        if (hs.find(p.first) != hs.end()) {
            counter = std::max(counter, hs[p.first]);
        }
    }
    return counter;
}
