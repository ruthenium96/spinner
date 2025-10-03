#include "Symmetrizer.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#ifndef NDEBUG
#include <stdexcept>
#include <string>
#endif
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace {
#ifndef NDEBUG
std::vector<double> squares_of_coefficients(
    size_t totalSpaceSize, 
    const std::vector<space::Subspace>& subspace) {
    std::vector<double> squares_of_coefficients(totalSpaceSize, 0);
    for (const auto& subspace : subspace) {
        for (size_t col = 0; col < subspace.decomposition->size_cols(); ++col) {
            double norm = 0.0;
            auto iterator = subspace.decomposition->GetNewIterator(col);
            while (iterator->hasNext()) {
                auto pair = iterator->getNext();
                norm += pair.value * pair.value;
            }
            iterator = subspace.decomposition->GetNewIterator(col);
            while (iterator->hasNext()) {
                auto pair = iterator->getNext();
                squares_of_coefficients[pair.index] += pair.value * pair.value / norm;
            }
        }
    }
    return squares_of_coefficients;
}

void compare_two_vectors_of_squared_coefficients(
    const std::vector<double>& before,
    const std::vector<double>& after) {
    std::vector<std::pair<uint32_t, double>> non_zero_differences;
    for (int i = 0; i < after.size(); ++i) {
        double value_diff = after[i] - before[i];
        if (std::abs(value_diff) > 1e-9) {
            non_zero_differences.push_back({i, value_diff});
        }
    }
    if (!non_zero_differences.empty()) {
        std::string error_message = "Symmetrizer: transformation is not unitary, indexes of problem and value_diff:\n";
        for (const auto& pair : non_zero_differences) {
            error_message += std::to_string(pair.first) + " : " + std::to_string(pair.second) + "\n";
        }
        throw std::logic_error(error_message);
    }
}
#endif

} // namespace

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
#ifndef NDEBUG
    std::vector<double> squares_of_coefficients_before = squares_of_coefficients(totalSpaceSize, space.getBlocks());
#endif
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

        // Auxiliary hash tables. They help to do not check orthogonality over and over.
        std::vector<std::unordered_map<uint32_t, std::vector<size_t>>> indexes_to_vectors_maps(
            group_.properties.number_of_representations);

        for (uint32_t l = 0; l < subspace_parent.decomposition->size_cols(); ++l) {
            uint8_t dimension_of_parent = subspace_parent.properties.dimensionality;
            // when we work with basi, we actually do all work for the orbit of basi,
            // so we add basi and its orbits to visited,
            // because there is no reason to work with them over and over
            std::vector<
                std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>>
                projected_basi = get_symmetrical_projected_decompositions(subspace_parent, l);
            for (size_t repr = 0; repr < group_.properties.number_of_representations; ++repr) {
                size_t j = group_.properties.number_of_representations * i + repr;
                gram_schmidt_orthogonalize(
                    projected_basi[repr], 
                    indexes_to_vectors_maps[repr], 
                    vector_result[j]);
                for (size_t k = 0;
                    k < group_.properties.number_of_projectors_of_representation[repr];
                    ++k) {
                    if (projected_basi[repr]->vempty(k)) {
                        // check if the DecompositionMap is empty:
                        continue;
                    }
                    move_vector_and_remember_it(
                        projected_basi[repr],
                        k,
                        indexes_to_vectors_maps[repr],
                        vector_result[j]);
                }
            }
        }
        subspace_parent.decomposition->clear();
    }

#ifndef NDEBUG
    std::vector<double> squares_of_coefficients_after = squares_of_coefficients(totalSpaceSize, vector_result);
    compare_two_vectors_of_squared_coefficients(squares_of_coefficients_before, squares_of_coefficients_after);
#endif
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

    for (auto& v : projections) {
        gram_schmidt_selforthogonalize(v);
        // this function deletes all pairs (key, value) with value = 0.0 from projections
        v->eraseExplicitZeros();
    }
    return projections;
}

void Symmetrizer::gram_schmidt_orthogonalize(
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>& decomposition_from,
    const std::unordered_map<uint32_t, std::vector<size_t>>& indexes_to_vectors_map,
    const Subspace& subspace_to) {
    std::unordered_set<uint32_t> indexes_in_decomposition;
    for (size_t col_i = 0; col_i < decomposition_from->size_cols(); ++col_i) {
        auto outer_iterator = decomposition_from->GetNewIterator(col_i);
        while (outer_iterator->hasNext()) {
            auto item = outer_iterator->getNext();
            indexes_in_decomposition.emplace(item.index);
        }
    }
    std::unordered_set<uint32_t> visited_vectors;
    for (const auto& index : indexes_in_decomposition) {
        if (!indexes_to_vectors_map.contains(index)) {
            continue;
        }
        for (const auto& col_j : indexes_to_vectors_map.at(index)) {
            if (visited_vectors.contains(col_j)) {
                continue;
            }
            auto first_iterator = subspace_to.decomposition->GetNewIterator(col_j);
            double norm = 0.0;
            while (first_iterator->hasNext()) {
                auto pair = first_iterator->getNext();
                norm += pair.value * pair.value;
            }
            for (size_t col_i = 0; col_i < decomposition_from->size_cols(); ++col_i) {
                auto second_iterator = decomposition_from->GetNewIterator(col_i);
                double dot_product = 0.0;
                while (second_iterator->hasNext()) {
                    auto pair = second_iterator->getNext();
                    auto index = pair.index;
                    auto value = pair.value;
                    if (!subspace_to.decomposition->is_zero(col_j, index)) {
                        dot_product += value * subspace_to.decomposition->at(col_j, index);
                    }
                }
                if (std::abs(dot_product) > 1e-9) {
                    auto third_iterator = subspace_to.decomposition->GetNewIterator(col_j);
                    while (third_iterator->hasNext()) {
                        auto pair = third_iterator->getNext();
                        double value = pair.value * dot_product / norm;
                        decomposition_from->add_to_position(-value, col_i, pair.index);
                    }
                }
            }
            visited_vectors.emplace(col_j);   
        }
    }
    decomposition_from->eraseExplicitZeros();
}

void Symmetrizer::gram_schmidt_selforthogonalize(
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&
        decomposition_from) {
    for (size_t col_j = 0; col_j < decomposition_from->size_cols(); ++col_j) {
        auto first_iterator = decomposition_from->GetNewIterator(col_j);
        double norm = 0.0;
        while (first_iterator->hasNext()) {
            auto pair = first_iterator->getNext();
            norm += pair.value * pair.value;
        }
        for (size_t col_i = col_j + 1; col_i < decomposition_from->size_cols(); ++col_i) {
            auto second_iterator = decomposition_from->GetNewIterator(col_i);
            double dot_product = 0.0;
            while (second_iterator->hasNext()) {
                auto pair = second_iterator->getNext();
                auto index = pair.index;
                auto value = pair.value;
                if (!decomposition_from->is_zero(col_j, index)) {
                    dot_product += value * decomposition_from->at(col_j, index);
                }
            }
            if (std::abs(dot_product) > 1e-9) {
                auto third_iterator = decomposition_from->GetNewIterator(col_j);
                while (third_iterator->hasNext()) {
                    auto pair = third_iterator->getNext();
                    double value = pair.value * dot_product / norm;
                    decomposition_from->add_to_position(-value, col_i, pair.index);
                }
            }
        }
    }
    decomposition_from->eraseExplicitZeros();
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