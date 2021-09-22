#include "Symmetrizer.h"

#include <utility>

Symmetrizer::Symmetrizer(spaces::LexicographicIndexConverter converter, Group group)
: converter_(std::move(converter)), group_(std::move(group)) {
}

Space Symmetrizer::apply(Space& space) const {
    std::vector<Subspace> vector_result;
    vector_result.resize(space.blocks.size() * group_.info.number_of_representations);
    entities::Entity::History history_result = space.history;

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.blocks.size(); ++i) {
        Subspace& subspace_parent = space.blocks[i];
        // add child subspaces of all representation (even if they will be empty)
        for (size_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
            BlockProperties block_properties = subspace_parent.properties;
            block_properties.representation.emplace_back(repr);
            block_properties.dimensionality *= group_.info.dimension_of_representation[repr];
            vector_result[group_.info.number_of_representations * i + repr].properties = block_properties;
        }

        // It is an auxiliary hash table. It helps to calculate each orbit only "dimensionality" times (see below).
        std::unordered_map<uint32_t, uint8_t> visited;
        // It is also an auxiliary hash table. It helps to do not check orthogonality over and over.
        std::vector<std::unordered_map<uint32_t, std::vector<size_t>>> added (group_.info.number_of_representations);

        for (auto & basi : subspace_parent.basis) {
            uint8_t dimension_of_parent = subspace_parent.properties.dimensionality;
            // when we work with basi, we actually do all work for the orbit of basi,
            // so we add basi and its orbits to visited,
            // because there is no reason to work with them over and over
            if (count_in_hash_table(basi, visited) < dimension_of_parent) {
                std::vector<std::vector<DecompositionMap>> projected_basi = get_symmetrical_projected_decompositions(basi);
                increment_in_hash_table(projected_basi[0][0], visited);
                for (size_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
                    for (size_t k = 0; k < group_.info.number_of_projectors_of_representation[repr]; ++k) {
                        if (projected_basi[repr][k].empty()) {
                            // check if the DecompositionMap is empty:
                            continue;
                        }
                        size_t j = group_.info.number_of_representations * i + repr;
                        add_vector_if_orthogonal_to_others(projected_basi[repr][k], added[repr], vector_result[j].basis);
                    }
                }
            }
        }

        subspace_parent.basis.clear();
    }

    return Space(std::move(vector_result), history_result);
}

std::vector<std::vector<DecompositionMap>> Symmetrizer::get_symmetrical_projected_decompositions(DecompositionMap & m) const {
    // it is a set (partitioned by representations) of all projected decompositions:
    std::vector<std::vector<DecompositionMap>> projections(group_.info.number_of_representations);
    for (uint8_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
        projections[repr].resize(group_.info.number_of_projectors_of_representation[repr]);
    }

    for (auto& p : m) {
        std::vector<uint8_t> nzs = converter_.convert_lex_index_to_sz_projections(p.first);
        std::vector<std::vector<uint8_t>> permutated_vectors = group_.permutate(nzs);

        for (uint8_t g = 0; g < group_.info.group_size; ++g) {
            uint32_t permutated_lex = converter_.convert_sz_projections_to_lex_index(permutated_vectors[g]);
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

void Symmetrizer::increment_in_hash_table(DecompositionMap & m,
                                          std::unordered_map<uint32_t , uint8_t>& hs) {
    for (const auto p : m) {
        if (hs.find(p.first) == hs.end()) {
            hs[p.first] = 1;
        } else {
            ++hs[p.first];
        }
    }
}

void Symmetrizer::erase_if_zero(std::vector<std::vector<DecompositionMap>>& projections) {
    for (auto & v : projections) {
        for (DecompositionMap& mm : v) {
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

uint8_t Symmetrizer::count_in_hash_table(const DecompositionMap & m,
                                   std::unordered_map<uint32_t , uint8_t>& hs) {
    uint8_t counter = 0;
    for (const auto p : m) {
        if (hs.find(p.first) != hs.end()) {
            counter = std::max(counter, hs[p.first]);
        }
    }
    return counter;
}

void Symmetrizer::add_vector_if_orthogonal_to_others(DecompositionMap &m,
                                                     std::unordered_map<uint32_t, std::vector<size_t>> &hs,
                                                     std::vector<DecompositionMap> &basis) {
    // TODO: should we check this only once per orbit?
    std::unordered_set<size_t> us;
    // we want to check orthogonality only with vectors, including the same lex-vectors:
    for (const auto p : m) {
        // hs[p.first] -- all vectors, including p.first lex-vector:
        for (const auto& lex : hs[p.first]) {
            // we do not want to check vector twice (or more):
            if (us.count(lex) > 0) {
                continue;
            }
            double accumulator = 0;
            for (const auto pp : m) {
                if (basis[lex].find(pp.first) != basis[lex].end()) {
                    accumulator += pp.second * basis[lex][pp.first];
                }
            }
            if (accumulator != 0) {
                return;
            }
            us.insert(lex);
        }
    }
    // if we reach this line -- DecompositionMap is okay, we can add it
    basis.emplace_back(std::move(m));
    for (auto p : basis.back()) {
        hs[p.first].emplace_back(basis.size() - 1);
    }
}
