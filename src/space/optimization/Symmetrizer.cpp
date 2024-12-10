#include "Symmetrizer.h"

#include <cassert>
#include <utility>

namespace space::optimization {

Symmetrizer::Symmetrizer(
    std::shared_ptr<const index_converter::AbstractIndexPermutator> permutator,
    group::Group group,
    quantum::linear_algebra::FactoriesList factories) :
    permutator_(std::move(permutator)),
    group_(std::move(group)),
    factories_(std::move(factories)) {}

Space Symmetrizer::apply(Space&& space) const {
    assert(!space.getBlocks().empty());
    auto totalSpaceSize = space.getBlocks()[0].decomposition->size_rows();
    std::vector<Subspace> vector_result;
    vector_result.reserve(space.getBlocks().size() * group_.properties.number_of_representations);
    for (size_t i = 0; i < space.getBlocks().size() * group_.properties.number_of_representations;
         ++i) {
        vector_result.emplace_back(factories_.createSparseSemiunitaryMatrix(0, totalSpaceSize));
    }

#pragma omp parallel for shared(space, vector_result) default(none)
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace_parent = space.getBlocks()[i];
        // add child subspaces of all representation (even if they will be empty)
        for (size_t repr = 0; repr < group_.properties.number_of_representations; ++repr) {
            BlockProperties block_properties = subspace_parent.properties;
            block_properties.representation.emplace_back(repr);
            block_properties.dimensionality *= group_.properties.dimension_of_representation[repr];
            vector_result[group_.properties.number_of_representations * i + repr].properties =
                block_properties;
        }

        // It is an auxiliary hash table. It helps to calculate each orbit only "dimensionality" times (see below).
        std::unordered_map<uint32_t, uint8_t> visited;
        // It is also an auxiliary hash table. It helps to do not check orthogonality over and over.
        std::vector<std::unordered_map<uint32_t, std::vector<size_t>>> added(
            group_.properties.number_of_representations);

        for (uint32_t l = 0; l < subspace_parent.decomposition->size_cols(); ++l) {
            uint8_t dimension_of_parent = subspace_parent.properties.dimensionality;
            // when we work with basi, we actually do all work for the orbit of basi,
            // so we add basi and its orbits to visited,
            // because there is no reason to work with them over and over
            if (count_how_many_orbit_was_visited(subspace_parent.decomposition, l, visited)
                < dimension_of_parent) {
                std::vector<
                    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>>
                    projected_basi = get_symmetrical_projected_decompositions(subspace_parent, l);
                increment_visited(projected_basi[0], 0, visited);
                for (size_t repr = 0; repr < group_.properties.number_of_representations; ++repr) {
                    for (size_t k = 0;
                         k < group_.properties.number_of_projectors_of_representation[repr];
                         ++k) {
                        if (projected_basi[repr]->vempty(k)) {
                            // check if the DecompositionMap is empty:
                            continue;
                        }
                        size_t j = group_.properties.number_of_representations * i + repr;
                        if (is_orthogonal_to_others(
                                projected_basi[repr],
                                k,
                                added[repr],
                                vector_result[j])) {
                            move_vector_and_remember_it(
                                projected_basi[repr],
                                k,
                                added[repr],
                                vector_result[j]);
                        }
                    }
                }
            }
        }

        subspace_parent.decomposition->clear();
    }

    return Space(std::move(vector_result));
}

std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>>
Symmetrizer::get_symmetrical_projected_decompositions(Subspace& subspace, uint32_t index_of_vector)
    const {
    // it is a set (partitioned by representations) of all projected decompositions:
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>>
        projections;
    for (uint8_t repr = 0; repr < group_.properties.number_of_representations; ++repr) {
        projections.emplace_back(std::move(factories_.createSparseSemiunitaryMatrix(
            group_.properties.number_of_projectors_of_representation[repr],
            permutator_->get_total_space_size())));
    }

    auto iterator = subspace.decomposition->GetNewIterator(index_of_vector);
    while (iterator->hasNext()) {
        auto item = iterator->getNext();
        auto permutated_indexes_and_signs = 
            permutator_->convert_index_to_permutated_indexes(item.index);

        for (uint8_t g = 0; g < group_.properties.group_size; ++g) {
            uint32_t permutated_index = permutated_indexes_and_signs[g].index;
            for (uint8_t repr = 0; repr < group_.properties.number_of_representations; ++repr) {
                for (uint8_t projector = 0;
                     projector < group_.properties.number_of_projectors_of_representation[repr];
                     ++projector) {
                    double value = group_.properties.coefficients_of_projectors[repr][projector][g]
                        * item.value * permutated_indexes_and_signs[g].sign;
                    projections[repr]->add_to_position(
                        value,
                        projector,
                        permutated_index);
                }
            }
        }
    }

    // this function deletes all pairs (key, value) with value = 0.0 from projections
    for (auto& v : projections) {
        v->eraseExplicitZeros();
    }
    return projections;
}

void Symmetrizer::increment_visited(
    const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& decomposition,
    uint32_t index_of_vector,
    std::unordered_map<uint32_t, uint8_t>& hs) {
    auto iterator = decomposition->GetNewIterator(index_of_vector);
    while (iterator->hasNext()) {
        auto item = iterator->getNext();
        if (hs.find(item.index) == hs.end()) {
            hs[item.index] = 1;
        } else {
            ++hs[item.index];
        }
    }
}

uint8_t Symmetrizer::count_how_many_orbit_was_visited(
    const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& decomposition,
    uint32_t index_of_vector,
    std::unordered_map<uint32_t, uint8_t>& hs) {
    uint8_t maximum = 0;
    auto iterator = decomposition->GetNewIterator(index_of_vector);
    while (iterator->hasNext()) {
        auto item = iterator->getNext();
        if (hs.find(item.index) != hs.end()) {
            maximum = std::max(maximum, hs[item.index]);
        }
    }
    return maximum;
}

bool Symmetrizer::is_orthogonal_to_others(
    const std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
        decomposition_from,
    uint32_t index_of_vector,
    std::unordered_map<uint32_t, std::vector<size_t>>& hs,
    const Subspace& subspace_to) {
    // TODO: should we check this only once per orbit?
    std::unordered_set<size_t> us;
    // we want to check orthogonality only with vectors, including the same lex-vectors:
    auto outer_iterator = decomposition_from->GetNewIterator(index_of_vector);
    while (outer_iterator->hasNext()) {
        auto item = outer_iterator->getNext();
        uint32_t index = item.index;
        double value = item.value;
        for (const auto& lex : hs[index]) {
            // we do not want to check vector twice (or more):
            if (us.count(lex) > 0) {
                continue;
            }
            double accumulator = 0;
            auto inner_iterator = decomposition_from->GetNewIterator(index_of_vector);
            while (inner_iterator->hasNext()) {
                auto inner_item = inner_iterator->getNext();
                uint32_t inner_index = inner_item.index;
                double inner_value = inner_item.value;
                if (!subspace_to.decomposition->is_zero(lex, inner_index)) {
                    accumulator += inner_value * subspace_to.decomposition->at(lex, inner_index);
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

void Symmetrizer::move_vector_and_remember_it(
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& decomposition_from,
    uint32_t index_of_vector,
    std::unordered_map<uint32_t, std::vector<size_t>>& hs,
    Subspace& subspace_to) {
    // if we reach this line -- DecompositionMap is okay, we can add it
    subspace_to.decomposition->move_vector_from(index_of_vector, decomposition_from);
    auto iterator =
        subspace_to.decomposition->GetNewIterator(subspace_to.decomposition->size_cols() - 1);
    while (iterator->hasNext()) {
        auto item = iterator->getNext();
        hs[item.index].emplace_back(subspace_to.decomposition->size_cols() - 1);
    }
}
}  // namespace space::optimization