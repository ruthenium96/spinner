#include "FTLMEigendecompositor.h"

#include <cassert>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "src/common/Quantity.h"
#include "src/eigendecompositor/ExactEigendecompositor.h"
#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"
#include "src/entities/spectrum/Subspectrum.h"

namespace eigendecompositor {

FTLMEigendecompositor::FTLMEigendecompositor(
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
    quantum::linear_algebra::FactoriesList factories_list,
    size_t krylov_subspace_size,
    size_t exact_decomposition_threshold,
    size_t number_of_seeds) :
    ExactEigendecompositor(converter, factories_list),
    converter_(converter),
    factories_list_(std::move(factories_list)),
    krylov_subspace_size_(krylov_subspace_size),
    exact_decomposition_threshold_(exact_decomposition_threshold),
    number_of_seeds_(number_of_seeds) {}
    
std::optional<OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
FTLMEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {

    uint32_t size_of_subspace = subspace.size();
    if (size_of_subspace <= exact_decomposition_threshold_) {
        if (!first_iteration_has_been_done_) {
            auto vector = factories_list_.createVector();
            vector->add_identical_values(size_of_subspace, 1.0);

            Subspectrum squared_back_projection_subspectrum;
            squared_back_projection_subspectrum.properties = subspace.properties;
            squared_back_projection_subspectrum.raw_data = std::move(vector);

            squared_back_projection_spectrum_[number_of_block] = std::move(squared_back_projection_subspectrum);
        }
        return ExactEigendecompositor::BuildSubspectra(number_of_block, subspace);
    }

    if (!first_iteration_has_been_done_) {
        seed_vectors_[number_of_block] = factories_list_.createRandomUnitVectors(size_of_subspace, number_of_seeds_);
        for (int seed = 0; seed < number_of_seeds_; ++seed) {
            weights_[number_of_block][seed] = factories_list_.createVector();
            weights_[number_of_block][seed]->add_identical_values(krylov_subspace_size_, subspace.properties.degeneracy * size_of_subspace);        
        }
    }

    std::optional<std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
        mb_unitary_transformation_matrix;

    // return_sparse_if_possible is true, because krylov eigendecomposition of sparse matrix is faster
    auto hamiltonian_submatrix = Submatrix(subspace, *energy_operator_, converter_, factories_list_, true);

    if (!do_we_need_eigenvectors_) {
        std::vector<Subspectrum> squared_back_projection_subspectrums;
        squared_back_projection_subspectrums.resize(number_of_seeds_);
        // if we need to explicitly calculate _only_ energy, we do not need eigenvectors:
        for (int seed = 0; seed < number_of_seeds_; ++seed) {
            auto krylov_couple = 
                hamiltonian_submatrix.raw_data->krylovDiagonalizeValues(
                    seed_vectors_[number_of_block][seed], 
                    krylov_subspace_size_);
        
            Subspectrum energy_subspectrum, squared_back_projection_subspectrum;
            energy_subspectrum.raw_data = std::move(krylov_couple.eigenvalues);
            energy_subspectrum.properties = hamiltonian_submatrix.properties;
            energy_subspectrum.properties.degeneracy *= (double)size_of_subspace;
    
            squared_back_projection_subspectrum.raw_data = std::move(krylov_couple.squared_back_projection);
            squared_back_projection_subspectrum.properties = hamiltonian_submatrix.properties;
            squared_back_projection_subspectrum.properties.degeneracy *= (double)size_of_subspace;
    
            energy_spectra_[number_of_block][seed] = std::move(energy_subspectrum);
            squared_back_projection_subspectrums[seed] = std::move(squared_back_projection_subspectrum);
        }
        squared_back_projection_spectrum_[number_of_block] = std::move(squared_back_projection_subspectrums);
    } else {
        mb_unitary_transformation_matrix = 
            std::vector<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>(number_of_seeds_);
        
        std::vector<Subspectrum> squared_back_projection_subspectrums;
        squared_back_projection_subspectrums.resize(number_of_seeds_);
        for (int seed = 0; seed < number_of_seeds_; ++seed) {
            auto krylov_triple = 
                hamiltonian_submatrix.raw_data->krylovDiagonalizeValuesVectors(
                    seed_vectors_[number_of_block][seed], 
                    krylov_subspace_size_);
        
            Subspectrum energy_subspectrum, squared_back_projection_subspectrum;
            energy_subspectrum.raw_data = std::move(krylov_triple.eigenvalues);
            energy_subspectrum.properties = hamiltonian_submatrix.properties;
            energy_subspectrum.properties.degeneracy *= (double)size_of_subspace;
    
            squared_back_projection_subspectrum.raw_data = std::move(krylov_triple.squared_back_projection);
            squared_back_projection_subspectrum.properties = hamiltonian_submatrix.properties;
            squared_back_projection_subspectrum.properties.degeneracy *= (double)size_of_subspace;
    
            energy_spectra_[number_of_block][seed] = std::move(energy_subspectrum);
            squared_back_projection_subspectrums[seed] = std::move(squared_back_projection_subspectrum);
            mb_unitary_transformation_matrix.value()[seed] = std::move(krylov_triple.eigenvectors);
        }
        squared_back_projection_spectrum_[number_of_block] = std::move(squared_back_projection_subspectrums);
    }
#ifndef NDEBUG
    energy_matrix_[number_of_block] = std::move(hamiltonian_submatrix);
#endif
    return mb_unitary_transformation_matrix;
}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
FTLMEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == common::Energy) {
        auto exact_subspectrumref = ExactEigendecompositor::getSubspectrum(common::Energy, number_of_block).value();

        const auto& ftlm_subspectrumref = energy_spectra_[number_of_block];

        if (getOneRef(exact_subspectrumref).get().raw_data != nullptr) {
            if (ftlm_subspectrumref[0].raw_data != nullptr) {
                throw std::logic_error("Both exact and FTLM constructed Energy block #" + std::to_string(number_of_block));
            }
            return exact_subspectrumref;
        } else {
            if (ftlm_subspectrumref[0].raw_data == nullptr) {
                throw std::logic_error("Neither exact nor FTLM constructed Energy block #" + std::to_string(number_of_block));
            }
            std::vector<std::reference_wrapper<const Subspectrum>> answer;
            for (const auto& el : ftlm_subspectrumref) {
                answer.push_back(std::reference_wrapper(el));
            }
            return answer;
        }
    }
    return std::nullopt;
}

std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
FTLMEigendecompositor::getSubmatrix(common::QuantityEnum quantity_enum, size_t number_of_block) const {
#ifndef NDEBUG
    if (quantity_enum == common::Energy) {
        auto exact_submatrixref = ExactEigendecompositor::getSubmatrix(common::Energy, number_of_block).value();

        const auto& ftlm_submatrixref = energy_matrix_[number_of_block];

        if (getOneRef(exact_submatrixref).get().raw_data != nullptr) {
            if (ftlm_submatrixref.raw_data != nullptr) {
                throw std::logic_error("Both exact and FTLM constructed Energy block #" + std::to_string(number_of_block));
            }
            return exact_submatrixref;
        } else {
            if (ftlm_submatrixref.raw_data == nullptr) {
                throw std::logic_error("Neither exact nor FTLM constructed Energy block #" + std::to_string(number_of_block));
            }
            return ftlm_submatrixref;
        }
    }
#endif
    return std::nullopt;
}

OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>
FTLMEigendecompositor::getWeightsOfBlockStates(size_t number_of_block) const {
    auto exact_weights = ExactEigendecompositor::getWeightsOfBlockStates(number_of_block);

    const auto& ftlm_weights = weights_[number_of_block];

    if (getOneRef(exact_weights).get() != nullptr) {
        if (ftlm_weights[0] != nullptr) {
            throw std::logic_error("Both exact and FTLM constructed weights of block #" + std::to_string(number_of_block));
        }
        return exact_weights;
    } else {
        if (ftlm_weights[0] == nullptr) {
            throw std::logic_error("Neither exact nor FTLM constructed weights of block #" + std::to_string(number_of_block));
        }
        std::vector<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> answer;
        for (const auto& el : ftlm_weights) {
            answer.push_back(std::reference_wrapper(el));
        }
        return answer;
    }
}

void FTLMEigendecompositor::initialize(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    uint32_t number_of_subspaces) {
    energy_spectra_.clear();
    energy_spectra_.resize(number_of_subspaces);
    for (int i = 0; i < number_of_subspaces; ++i) {
        energy_spectra_[i].clear();
        energy_spectra_[i].resize(number_of_seeds_);
    }

    energy_matrix_.clear();
    energy_matrix_.resize(number_of_subspaces);

    if (!first_iteration_has_been_done_) {
        seed_vectors_.resize(number_of_subspaces);
        squared_back_projection_spectrum_.resize(number_of_subspaces);
        weights_.resize(number_of_subspaces);
        for (int i = 0; i < number_of_subspaces; ++i) {
            weights_[i].clear();
            weights_[i].resize(number_of_seeds_);
        }
    }

    energy_operator_ = operators_to_calculate.at(common::Energy);
    do_we_need_eigenvectors_ =
        !(operators_to_calculate.size() == 1 && derivatives_operators_to_calculate.empty());
    // Do not do
    // std::erase_if(operators_to_calculate, [](const auto& p) { return p.first == common::Energy; });
    // because ExactEigendecompositor need Hamiltonian too.
    ExactEigendecompositor::initialize(operators_to_calculate, derivatives_operators_to_calculate, number_of_subspaces);
}

void FTLMEigendecompositor::finalize() {
    ExactEigendecompositor::finalize();
    first_iteration_has_been_done_ = true;
}

}