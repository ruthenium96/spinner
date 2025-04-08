#include "FTLMEigendecompositor.h"

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include "src/common/Quantity.h"
#include "src/eigendecompositor/ExactEigendecompositor.h"
#include "src/entities/data_structures/AbstractDiagonalizableMatrix.h"
#include "src/entities/spectrum/Spectrum.h"
#include "src/entities/spectrum/Subspectrum.h"

namespace eigendecompositor {

FTLMEigendecompositor::FTLMEigendecompositor(
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
    quantum::linear_algebra::FactoriesList factories_list,
    size_t krylov_subspace_size,
    size_t exact_decomposition_threshold) :
    ExactEigendecompositor(converter, factories_list),
    converter_(converter),
    factories_list_(std::move(factories_list)),
    krylov_subspace_size_(krylov_subspace_size),
    exact_decomposition_threshold_(exact_decomposition_threshold) {}
    
std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
FTLMEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {

    uint32_t size_of_subspace = subspace.size();
    if (size_of_subspace <= exact_decomposition_threshold_) {
        if (squared_back_projection_.spectrum_.blocks[number_of_block].raw_data == nullptr) {
            auto vector = factories_list_.createVector();
            vector->add_identical_values(size_of_subspace, 1.0);

            Subspectrum squared_back_projection_subspectrum;
            squared_back_projection_subspectrum.properties = subspace.properties;
            squared_back_projection_subspectrum.raw_data = std::move(vector);

            squared_back_projection_.spectrum_.blocks[number_of_block] = std::move(squared_back_projection_subspectrum);
        }
        return ExactEigendecompositor::BuildSubspectra(number_of_block, subspace);
    }

    if (seed_vectors_.at(number_of_block) == nullptr) {
        seed_vectors_[number_of_block] = std::move(factories_list_.createRandomUnitVector(size_of_subspace));
    }

    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        mb_unitary_transformation_matrix;

    // return_sparse_if_possible is true, because krylov eigendecomposition of sparse matrix is faster
    auto hamiltonian_submatrix = Submatrix(subspace, *energy_operator_, converter_, factories_list_, true);

    if (!do_we_need_eigenvectors_) {
        // if we need to explicitly calculate _only_ energy, we do not need eigenvectors:
        auto krylov_couple = 
            hamiltonian_submatrix.raw_data->krylovDiagonalizeValues(seed_vectors_[number_of_block], krylov_subspace_size_);
        
        Subspectrum energy_subspectrum, squared_back_projection_subspectrum;
        energy_subspectrum.raw_data = std::move(krylov_couple.eigenvalues);
        energy_subspectrum.properties = hamiltonian_submatrix.properties;
        energy_subspectrum.properties.degeneracy *= (double)size_of_subspace;

        squared_back_projection_subspectrum.raw_data = std::move(krylov_couple.squared_back_projection);
        squared_back_projection_subspectrum.properties = hamiltonian_submatrix.properties;
        squared_back_projection_subspectrum.properties.degeneracy *= (double)size_of_subspace;

        energy_.spectrum_.blocks[number_of_block] = std::move(energy_subspectrum);
        squared_back_projection_.spectrum_.blocks[number_of_block] = std::move(squared_back_projection_subspectrum);
    } else {
        throw std::invalid_argument("NOT IMPLEMENTED YET!");
    }
#ifndef NDEBUG
    energy_.matrix_.blocks[number_of_block] = std::move(hamiltonian_submatrix);
#endif

    return mb_unitary_transformation_matrix;
}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
FTLMEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == common::Energy) {
        auto exact_subspectrumref = ExactEigendecompositor::getSubspectrum(common::Energy, number_of_block).value();

        const auto& ftlm_subspectrumref = energy_.spectrum_.blocks[number_of_block];

        if (std::get<std::reference_wrapper<const Subspectrum>>(exact_subspectrumref).get().raw_data != nullptr) {
            if (ftlm_subspectrumref.raw_data != nullptr) {
                throw std::logic_error("Both exact and FTLM constructed Energy block #" + std::to_string(number_of_block));
            }
            return exact_subspectrumref;
        } else {
            if (ftlm_subspectrumref.raw_data == nullptr) {
                throw std::logic_error("Neither exact or FTLM constructed Energy block #" + std::to_string(number_of_block));
            }
            return ftlm_subspectrumref;
        }
    }
    if (quantity_enum == common::squared_back_projection) {
        return squared_back_projection_.spectrum_.blocks[number_of_block];
    }
    return std::nullopt;
}

std::optional<MatrixRef>
FTLMEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
#ifndef NDEBUG
    if (quantity_enum == common::Energy) {
        MatrixRef exact_matrixref = ExactEigendecompositor::getMatrix(common::Energy).value();

        MatrixRef ftlm_matrixref = MatrixRef(energy_.matrix_);

        assert(underlying_matrixref.blocks.size() == ftlm_matrixref.blocks.size());
        for (int i = 0; i < ftlm_matrixref.blocks.size(); ++i) {
            if (exact_matrixref.blocks[i].get().raw_data != nullptr) {
                if (ftlm_matrixref.blocks[i].get().raw_data != nullptr) {
                    throw std::logic_error("Both exact and FTLM constructed Energy block");
                }
                ftlm_matrixref.blocks[i] = exact_matrixref.blocks[i];
            } else {
                if (ftlm_matrixref.blocks[i].get().raw_data == nullptr) {
                    throw std::logic_error("Neither exact or FTLM constructed Energy block");
                }
                // ftlm_matrixref.blocks[i] = ftlm_matrixref.blocks[i];
            }
        }
        return ftlm_matrixref;
    }
#endif
    return std::nullopt;
}

std::optional<SpectrumRef> FTLMEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return std::nullopt;
}

std::optional<MatrixRef> FTLMEigendecompositor::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return std::nullopt;
}

void FTLMEigendecompositor::initialize(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    uint32_t number_of_subspaces) {
    energy_.matrix_.blocks.clear();
    energy_.spectrum_.blocks.clear();
    squared_back_projection_.spectrum_.blocks.clear();
    energy_.matrix_.blocks.resize(number_of_subspaces);
    energy_.spectrum_.blocks.resize(number_of_subspaces);
    squared_back_projection_.spectrum_.blocks.resize(number_of_subspaces);

    if (seed_vectors_.empty()) {
        seed_vectors_.resize(number_of_subspaces);
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
}

}