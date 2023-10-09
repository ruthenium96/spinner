#include "ExplicitQuantitiesEigendecompositor.h"

#include <utility>

namespace eigendecompositor {

ExplicitQuantitiesEigendecompositor::ExplicitQuantitiesEigendecompositor(
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
    lexicographic::IndexConverter converter,
    quantum::linear_algebra::FactoriesList factories_list) :
    eigendecompositor_(std::move(eigendecompositor)),
    converter_(std::move(converter)),
    factories_list_(std::move(factories_list)) {}

std::optional<std::reference_wrapper<const Spectrum>>
ExplicitQuantitiesEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantities_map_.contains(quantity_enum)) {
        return quantities_map_.at(quantity_enum).spectrum_;
    }
    return eigendecompositor_->getSpectrum(quantity_enum);
}

std::optional<std::reference_wrapper<const Matrix>>
ExplicitQuantitiesEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantities_map_.contains(quantity_enum)) {
        return quantities_map_.at(quantity_enum).matrix_;
    }
    return eigendecompositor_->getMatrix(quantity_enum);
}

std::optional<std::reference_wrapper<const Spectrum>>
ExplicitQuantitiesEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (derivatives_map_.contains({quantity_enum, symbol_name})) {
        return derivatives_map_.at({quantity_enum, symbol_name}).spectrum_;
    }
    return eigendecompositor_->getSpectrumDerivative(quantity_enum, symbol_name);
}

std::optional<std::reference_wrapper<const Matrix>>
ExplicitQuantitiesEigendecompositor::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (derivatives_map_.contains({quantity_enum, symbol_name})) {
        return derivatives_map_.at({quantity_enum, symbol_name}).matrix_;
    }
    return eigendecompositor_->getMatrixDerivative(quantity_enum, symbol_name);
}

std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
ExplicitQuantitiesEigendecompositor::BuildSubspectra(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    size_t number_of_block,
    const space::Subspace& subspace) {
    auto mb_unitary_transformation_matrix = eigendecompositor_->BuildSubspectra(
        operators_to_calculate,
        derivatives_operators_to_calculate,
        number_of_block,
        subspace);

    for (const auto& [quantity_enum, operator_to_calculate] : operators_to_calculate) {
        if (quantity_enum == common::Energy) {
            throw std::invalid_argument(
                "Energy operator passed to ExplicitQuantitiesEigendecompositor");
        }
        auto non_hamiltonian_submatrix =
            Submatrix(subspace, *operator_to_calculate, converter_, factories_list_);
        if (!quantities_map_.contains(quantity_enum)) {
            quantities_map_[quantity_enum] = common::Quantity();
        }
        auto& quantity = quantities_map_[quantity_enum];
        quantity.spectrum_.blocks.emplace_back(non_energy_subspectrum(
            non_hamiltonian_submatrix,
            mb_unitary_transformation_matrix.value()));
        quantity.matrix_.blocks.emplace_back(std::move(non_hamiltonian_submatrix));
    }

    for (auto& [pair, derivative_operator] : derivatives_operators_to_calculate) {
        auto derivative_submatrix =
            Submatrix(subspace, *derivative_operator, converter_, factories_list_);
        if (!derivatives_map_.contains(pair)) {
            derivatives_map_[pair] = common::Quantity();
        }
        auto& derivative = derivatives_map_[pair];
        derivative.spectrum_.blocks.emplace_back(
            non_energy_subspectrum(derivative_submatrix, mb_unitary_transformation_matrix.value()));
        derivative.matrix_.blocks.emplace_back(std::move(derivative_submatrix));
    }

    // delete operators, because we just have calculated corresponding spectra
    for (const auto& p : quantities_map_) {
        auto quantity_enum = p.first;
        std::erase_if(operators_to_calculate, [quantity_enum](auto p) {
            return p.first == quantity_enum;
        });
    }
    for (const auto& p : derivatives_map_) {
        auto pair = p.first;
        std::erase_if(derivatives_operators_to_calculate, [pair](auto p) {
            return p.first == pair;
        });
    }

    return mb_unitary_transformation_matrix;
}

void ExplicitQuantitiesEigendecompositor::initialize() {
    for (auto& [_, quantity_] : quantities_map_) {
        quantity_.spectrum_.blocks.clear();
    }
    for (auto& [_, derivative] : derivatives_map_) {
        derivative.spectrum_.blocks.clear();
    }
    eigendecompositor_->initialize();
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