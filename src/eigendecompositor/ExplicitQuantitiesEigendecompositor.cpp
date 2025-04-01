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

std::optional<SpectrumRef>
ExplicitQuantitiesEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantities_map_.contains(quantity_enum)) {
        return SpectrumRef(quantities_map_.at(quantity_enum).spectrum_);
    }
    return eigendecompositor_->getSpectrum(quantity_enum);
}

std::optional<std::reference_wrapper<const Subspectrum>>
ExplicitQuantitiesEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantities_map_.contains(quantity_enum)) {
        return quantities_map_.at(quantity_enum).spectrum_.blocks[number_of_block];
    }
    return eigendecompositor_->getSubspectrum(quantity_enum, number_of_block);
}

std::optional<MatrixRef>
ExplicitQuantitiesEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantities_map_.contains(quantity_enum)) {
        return MatrixRef(quantities_map_.at(quantity_enum).matrix_);
    }
    return eigendecompositor_->getMatrix(quantity_enum);
}

std::optional<SpectrumRef>
ExplicitQuantitiesEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (derivatives_map_.contains({quantity_enum, symbol_name})) {
        return SpectrumRef(derivatives_map_.at({quantity_enum, symbol_name}).spectrum_);
    }
    return eigendecompositor_->getSpectrumDerivative(quantity_enum, symbol_name);
}

std::optional<MatrixRef>
ExplicitQuantitiesEigendecompositor::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (derivatives_map_.contains({quantity_enum, symbol_name})) {
        return MatrixRef(derivatives_map_.at({quantity_enum, symbol_name}).matrix_);
    }
    return eigendecompositor_->getMatrixDerivative(quantity_enum, symbol_name);
}

std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
ExplicitQuantitiesEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {
    auto mb_unitary_transformation_matrix =
        eigendecompositor_->BuildSubspectra(number_of_block, subspace);

    for (const auto& [quantity_enum, operator_to_calculate] : quantities_operators_map_) {
        auto non_hamiltonian_submatrix =
            Submatrix(subspace, *operator_to_calculate, converter_, factories_list_);
        auto& quantity = quantities_map_[quantity_enum];
        quantity.spectrum_.blocks[number_of_block] = non_energy_subspectrum(
            non_hamiltonian_submatrix,
            mb_unitary_transformation_matrix.value());
        quantity.matrix_.blocks[number_of_block] = std::move(non_hamiltonian_submatrix);
    }

    for (auto& [pair, derivative_operator] : derivatives_operators_map_) {
        auto derivative_submatrix =
            Submatrix(subspace, *derivative_operator, converter_, factories_list_);
        auto& derivative = derivatives_map_[pair];
        derivative.spectrum_.blocks[number_of_block] =
            non_energy_subspectrum(derivative_submatrix, mb_unitary_transformation_matrix.value());
        derivative.matrix_.blocks[number_of_block] = std::move(derivative_submatrix);
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
        if (!quantities_map_.contains(quantity_enum)) {
            quantities_map_[quantity_enum] = common::Quantity();
        }
        quantities_map_[quantity_enum].matrix_.blocks.clear();
        quantities_map_[quantity_enum].spectrum_.blocks.clear();
        quantities_map_[quantity_enum].matrix_.blocks.resize(number_of_subspaces);
        quantities_map_[quantity_enum].spectrum_.blocks.resize(number_of_subspaces);
    }

    for (auto& [pair, derivative_operator] : derivatives_operators_to_calculate) {
        derivatives_operators_map_[pair] = derivative_operator;
        if (!derivatives_map_.contains(pair)) {
            derivatives_map_[pair] = common::Quantity();
        }
        derivatives_map_[pair].matrix_.blocks.clear();
        derivatives_map_[pair].spectrum_.blocks.clear();
        derivatives_map_[pair].matrix_.blocks.resize(number_of_subspaces);
        derivatives_map_[pair].spectrum_.blocks.resize(number_of_subspaces);
    }

    // delete operators, because we are going to calculate corresponding spectra
    for (const auto& p : quantities_map_) {
        auto quantity_enum = p.first;
        std::erase_if(operators_to_calculate, [quantity_enum](const auto& p) {
            return p.first == quantity_enum;
        });
    }
    for (const auto& p : derivatives_map_) {
        auto pair = p.first;
        std::erase_if(derivatives_operators_to_calculate, [pair](const auto& p) {
            return p.first == pair;
        });
    }
}

void ExplicitQuantitiesEigendecompositor::finalize() {
    eigendecompositor_->finalize();
}

Subspectrum ExplicitQuantitiesEigendecompositor::non_energy_subspectrum(
    const Submatrix& non_hamiltonian_submatrix,
    const std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>&
        unitary_transformation_matrix) {
    auto raw_data = unitary_transformation_matrix->unitaryTransformAndReturnMainDiagonal(
        non_hamiltonian_submatrix.raw_data);

    auto non_energy_subspectrum =
        Subspectrum(std::move(raw_data), non_hamiltonian_submatrix.properties);

    return std::move(non_energy_subspectrum);
}
}  // namespace eigendecompositor