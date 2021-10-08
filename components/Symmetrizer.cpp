#include "Symmetrizer.h"

#include <utility>

namespace {
bool SizeOfPermutationsEqualsNumberOfSpins (const spaces::LexicographicIndexConverter& converter,
                                            const Group& group) {
    return converter.mults_.size() == group.elements_[0].size();
}

bool OrbitOfCentersHasTheSameValueOfMultiplicity (const spaces::LexicographicIndexConverter& converter,
                                                const Group& group) {
    for (const auto& el : group.elements_) {
        std::vector<int> permutated_mults(converter.mults_);
        for (size_t i = 0; i < group.elements_[0].size(); ++i) {
            permutated_mults[i] = converter.mults_[el[i]];
        }
        if (permutated_mults != converter.mults_) {
            return false;
        }
    }
    return true;
}
}

Symmetrizer::Symmetrizer(spaces::LexicographicIndexConverter converter, Group group)
: converter_(std::move(converter)), group_(std::move(group)) {
    if (!SizeOfPermutationsEqualsNumberOfSpins(converter_, group_)) {
        throw std::length_error("The size of group elements does not equal to the number of spins.");
    }
    if (!OrbitOfCentersHasTheSameValueOfMultiplicity(converter_, group_)) {
        throw std::invalid_argument("Group permutes centers with different multiplicities.");
    }
}

Space Symmetrizer::apply(Space& space) const {
    std::vector<Subspace> vector_result;
    vector_result.resize(space.blocks.size() * group_.info.number_of_representations);

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.blocks.size(); ++i) {
        Subspace& subspace_parent = space.blocks[i];
        // add child subspaces of all representation (even if they will be empty)
        for (size_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
            BlockProperties block_properties = subspace_parent.properties;
            block_properties.representation.emplace_back(repr);
            block_properties.dimensionality *= group_.info.dimension_of_representation[repr];
            vector_result[group_.info.number_of_representations * i + repr].properties = block_properties;
            // TODO: still bad idea
            vector_result[group_.info.number_of_representations * i + repr].tensor_size = subspace_parent.tensor_size;
        }

        // It is an auxiliary hash table. It helps to calculate each orbit only "dimensionality" times (see below).
        std::unordered_map<uint32_t, uint8_t> visited;
        // It is also an auxiliary hash table. It helps to do not check orthogonality over and over.
        std::vector<std::unordered_map<uint32_t, std::vector<size_t>>> added (group_.info.number_of_representations);

        for (uint32_t l = 0; l < subspace_parent.size(); ++l) {
            uint8_t dimension_of_parent = subspace_parent.properties.dimensionality;
            // when we work with basi, we actually do all work for the orbit of basi,
            // so we add basi and its orbits to visited,
            // because there is no reason to work with them over and over
            if (count_how_many_orbit_was_visited(subspace_parent, l, visited) < dimension_of_parent) {
                std::vector<Subspace> projected_basi = get_symmetrical_projected_decompositions(subspace_parent, l);
                increment_visited(projected_basi[0], 0, visited);
                for (size_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
                    for (size_t k = 0; k < group_.info.number_of_projectors_of_representation[repr]; ++k) {
                        if (projected_basi[repr].vempty(k)) {
                            // check if the DecompositionMap is empty:
                            continue;
                        }
                        size_t j = group_.info.number_of_representations * i + repr;
                        if (is_orthogonal_to_others(projected_basi[repr], k, added[repr], vector_result[j])) {
                            move_vector_and_remember_it(projected_basi[repr], k, added[repr], vector_result[j]);
                        }
                    }
                }
            }
        }

        subspace_parent.clear();
    }

    return Space(std::move(vector_result));
}

std::vector<Subspace> Symmetrizer::get_symmetrical_projected_decompositions(Subspace& subspace,
                                                                            uint32_t index_of_vector) const {
    // it is a set (partitioned by representations) of all projected decompositions:
    std::vector<Subspace> projections(group_.info.number_of_representations);
    for (uint8_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
        projections[repr].tensor_size = subspace.tensor_size;
        projections[repr].resize(group_.info.number_of_projectors_of_representation[repr]);
    }

    for (auto p = subspace.vbegin(index_of_vector); p != subspace.vend(index_of_vector); ++p) {
        std::vector<uint8_t> nzs = converter_.convert_lex_index_to_sz_projections(INDEX(p));
        std::vector<std::vector<uint8_t>> permutated_vectors = group_.permutate(nzs);

        for (uint8_t g = 0; g < group_.info.group_size; ++g) {
            uint32_t permutated_lex = converter_.convert_sz_projections_to_lex_index(permutated_vectors[g]);
            for (uint8_t repr = 0; repr < group_.info.number_of_representations; ++repr) {
                for (uint8_t projector = 0; projector < group_.info.number_of_projectors_of_representation[repr]; ++projector) {
                    projections[repr].add_to_position(group_.info.coefficients_of_projectors[repr][projector][g] * VALUE(p),
                                                      projector, permutated_lex);
                }
            }
        }
    }

    // this function deletes all pairs (key, value) with value = 0.0 from projections
    for (auto & v : projections) {
        v.erase_if_zero();
    }
    return projections;
}

void Symmetrizer::increment_visited(const Subspace& subspace,
                                    uint32_t index_of_vector,
                                    std::unordered_map<uint32_t , uint8_t>& hs) {
    for (auto p = subspace.vbegin(index_of_vector); p != subspace.vend(index_of_vector); ++p) {
        if (hs.find(INDEX(p)) == hs.end()) {
            hs[INDEX(p)] = 1;
        } else {
            ++hs[INDEX(p)];
        }
    }
}

uint8_t Symmetrizer::count_how_many_orbit_was_visited(const Subspace& subspace,
                                                      uint32_t index_of_vector,
                                                      std::unordered_map<uint32_t , uint8_t>& hs) {
    uint8_t maximum = 0;
    for (auto p = subspace.vbegin(index_of_vector); p != subspace.vend(index_of_vector); ++p) {
        if (hs.find(INDEX(p)) != hs.end()) {
            maximum = std::max(maximum, hs[INDEX(p)]);
        }

    }
    return maximum;
}

bool Symmetrizer::is_orthogonal_to_others(const Subspace& subspace_from, uint32_t index_of_vector,
                                          std::unordered_map<uint32_t, std::vector<size_t>>& hs,
                                          const Subspace& subspace_to) {
    // TODO: should we check this only once per orbit?
    std::unordered_set<size_t> us;
    // we want to check orthogonality only with vectors, including the same lex-vectors:
    for (auto p = subspace_from.vbegin(index_of_vector); p != subspace_from.vend(index_of_vector); ++p) {
        // hs[p.first] -- all vectors, including p.first lex-vector:
        for (const auto& lex : hs[INDEX(p)]) {
            // we do not want to check vector twice (or more):
            if (us.count(lex) > 0) {
                continue;
            }
            double accumulator = 0;
            for (auto pp = subspace_from.vbegin(index_of_vector); pp != subspace_from.vend(index_of_vector); ++pp) {
                if (!subspace_to.is_zero(lex, INDEX(pp))) {
                    accumulator += VALUE(pp) * subspace_to(lex, INDEX(pp));
                }
            }
            if (accumulator != 0) {
                return false;
            }
            us.insert(lex);
        }
    }
    return true;
}

void Symmetrizer::move_vector_and_remember_it(Subspace& subspace_from, uint32_t index_of_vector,
                                              std::unordered_map<uint32_t, std::vector<size_t>>& hs,
                                              Subspace& subspace_to) {
    // if we reach this line -- DecompositionMap is okay, we can add it
    subspace_to.move_vector_from(index_of_vector, subspace_from);
    for (auto p = subspace_to.vbegin(subspace_to.size() - 1); p != subspace_to.vend(subspace_to.size() - 1); ++p) {
        hs[INDEX(p)].emplace_back(subspace_to.size() - 1);
    }
}
