#include "ExplicitQuantitiesEigendecompositor.h"

#include <utility>

namespace eigendecompositor {

ExplicitQuantitiesEigendecompositor::ExplicitQuantitiesEigendecompositor(
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
    quantum::linear_algebra::FactoriesList factories_list) :
    eigendecompositor_(std::move(eigendecompositor)),
    converter_(std::move(converter)),
    factories_list_(std::move(factories_list)) {}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
ExplicitQuantitiesEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantities_spectra_map_.contains(quantity_enum)) {
        return copyRef<Subspectrum, std::reference_wrapper<const Subspectrum>>(
            quantities_spectra_map_.at(quantity_enum)[number_of_block]);
    }
    return eigendecompositor_->getSubspectrum(quantity_enum, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
ExplicitQuantitiesEigendecompositor::getSubmatrix(common::QuantityEnum quantity_enum, size_t number_of_block) const {
#ifndef NDEBUG
    if (quantities_matrix_map_.contains(quantity_enum)) {
        return quantities_matrix_map_.at(quantity_enum)[number_of_block];
    }
#endif
    return eigendecompositor_->getSubmatrix(quantity_enum, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
ExplicitQuantitiesEigendecompositor::getSubspectrumDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name, size_t number_of_block) const {
    if (derivatives_spectra_map_.contains({quantity_enum, symbol_name})) {
        return copyRef<Subspectrum, std::reference_wrapper<const Subspectrum>>(
            derivatives_spectra_map_.at({quantity_enum, symbol_name})[number_of_block]);
    }
    return eigendecompositor_->getSubspectrumDerivative(quantity_enum, symbol_name, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
ExplicitQuantitiesEigendecompositor::getSubmatrixDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name, size_t number_of_block) const {
#ifndef NDEBUG
    if (derivatives_matrix_map_.contains({quantity_enum, symbol_name})) {
        return derivatives_matrix_map_.at({quantity_enum, symbol_name})[number_of_block];
    }
#endif
    return eigendecompositor_->getSubmatrixDerivative(quantity_enum, symbol_name, number_of_block);
}

std::optional<OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
ExplicitQuantitiesEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {
    auto mb_unitary_transformation_matrix =
        eigendecompositor_->BuildSubspectra(number_of_block, subspace);

    for (const auto& [quantity_enum, operator_to_calculate] : quantities_operators_map_) {
        // return_sparse_if_possible is true, because unitary transformation of sparse matrix is faster
        auto non_hamiltonian_submatrix =
            Submatrix(subspace, *operator_to_calculate, converter_, factories_list_, true);
        auto& quantity_spectrum = quantities_spectra_map_[quantity_enum];
        quantity_spectrum[number_of_block] = non_energy_subspectrum(
            non_hamiltonian_submatrix,
            mb_unitary_transformation_matrix.value());
#ifndef NDEBUG
        auto& quantity_matrix = quantities_matrix_map_[quantity_enum];
        quantity_matrix[number_of_block] = std::move(non_hamiltonian_submatrix);
#endif
    }

    for (auto& [pair, derivative_operator] : derivatives_operators_map_) {
        // return_sparse_if_possible is true, because unitary transformation of sparse matrix is faster
        auto derivative_submatrix =
            Submatrix(subspace, *derivative_operator, converter_, factories_list_, true);
        auto& derivative_spectrum = derivatives_spectra_map_[pair];
        derivative_spectrum[number_of_block] =
            non_energy_subspectrum(derivative_submatrix, mb_unitary_transformation_matrix.value());
#ifndef NDEBUG
        auto& derivative_matrix = derivatives_matrix_map_[pair];
        derivative_matrix[number_of_block] = std::move(derivative_submatrix);
#endif
    }

    return mb_unitary_transformation_matrix;
}

void ExplicitQuantitiesEigendecompositor::initialize(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    uint32_t number_of_subspaces) {
    eigendecompositor_->initialize(
        operators_to_calculate,
        derivatives_operators_to_calculate,
        number_of_subspaces);

    for (const auto& [quantity_enum, operator_to_calculate] : operators_to_calculate) {
        if (quantity_enum == common::Energy) {
            throw std::invalid_argument(
                "Energy operator passed to ExplicitQuantitiesEigendecompositor");
        }
        quantities_operators_map_[quantity_enum] = operator_to_calculate;
        if (!quantities_spectra_map_.contains(quantity_enum)) {
            quantities_spectra_map_[quantity_enum] = std::vector<OneOrMany<Subspectrum>>();
        }
        if (!quantities_matrix_map_.contains(quantity_enum)) {
            quantities_matrix_map_[quantity_enum] = std::vector<Submatrix>();
        }
        quantities_matrix_map_[quantity_enum].clear();
        quantities_spectra_map_[quantity_enum].clear();
        quantities_matrix_map_[quantity_enum].resize(number_of_subspaces);
        quantities_spectra_map_[quantity_enum].resize(number_of_subspaces);
    }

    for (auto& [pair, derivative_operator] : derivatives_operators_to_calculate) {
        derivatives_operators_map_[pair] = derivative_operator;
        if (!derivatives_spectra_map_.contains(pair)) {
            derivatives_spectra_map_[pair] = std::vector<OneOrMany<Subspectrum>>();
        }
        if (!derivatives_matrix_map_.contains(pair)) {
            derivatives_matrix_map_[pair] = std::vector<Submatrix>();
        }
        derivatives_spectra_map_[pair].clear();
        derivatives_matrix_map_[pair].clear();
        derivatives_spectra_map_[pair].resize(number_of_subspaces);
        derivatives_matrix_map_[pair].resize(number_of_subspaces);
    }

    // delete operators, because we are going to calculate corresponding spectra
    for (const auto& p : quantities_spectra_map_) {
        auto quantity_enum = p.first;
        std::erase_if(operators_to_calculate, [quantity_enum](const auto& p) {
            return p.first == quantity_enum;
        });
    }
    for (const auto& p : derivatives_spectra_map_) {
        auto pair = p.first;
        std::erase_if(derivatives_operators_to_calculate, [pair](const auto& p) {
            return p.first == pair;
        });
    }
}

void ExplicitQuantitiesEigendecompositor::finalize() {
    eigendecompositor_->finalize();
}

OneOrMany<Subspectrum> ExplicitQuantitiesEigendecompositor::non_energy_subspectrum(
    const Submatrix& non_hamiltonian_submatrix,
    const OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>&
        unitary_transformation_matrix) {
    return transform_one_or_many(
    std::function([&non_hamiltonian_submatrix](std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix> unitary_transformation_matrix) {
        auto raw_data = unitary_transformation_matrix->unitaryTransformAndReturnMainDiagonal(
            non_hamiltonian_submatrix.raw_data);
    
        auto non_energy_subspectrum =
            Subspectrum(std::move(raw_data), non_hamiltonian_submatrix.properties);
    
        return std::move(non_energy_subspectrum);
    }), unitary_transformation_matrix);
}
}  // namespace eigendecompositor