#include "OneSymbolInHamiltonianEigendecompositor.h"

#include <cassert>
#include <utility>

namespace eigendecompositor {

OneSymbolInHamiltonianEigendecompositor::OneSymbolInHamiltonianEigendecompositor(
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
    std::function<double()> currentValueGetter) :
    eigendecompositor_(std::move(eigendecompositor)),
    currentValueGetter_(std::move(currentValueGetter)) {
    initial_value_of_symbol_ = currentValueGetter_();
    if (std::abs(initial_value_of_symbol_) < 1e-9) {
        throw std::invalid_argument(
            "Initial value of the symbol is too small to be used "
            "inside OneSymbolInHamiltonianEigendecompositor");
    }
}

std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
OneSymbolInHamiltonianEigendecompositor::BuildSubspectra(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    size_t number_of_block,
    const space::Subspace& subspace) {
    size_t number_of_all_derivatives = derivatives_operators_to_calculate.size();
    // delete energy derivative operator, because we are about to calculate corresponding spectrum
    std::erase_if(derivatives_operators_to_calculate, [](auto p) {
        return p.first.first == common::Energy;
    });
    size_t number_of_non_energy_derivatives = derivatives_operators_to_calculate.size();
    size_t number_of_energy_derivatives =
        number_of_all_derivatives - number_of_non_energy_derivatives;

    if (number_of_energy_derivatives >= 2) {
        throw std::invalid_argument(
            "More than two energy derivatives passed to OneParameterEigendecompositor");
    }
    if (number_of_energy_derivatives == 1 && !current_energy_derivative_spectrum_.has_value()) {
        current_energy_derivative_spectrum_ = Spectrum();
    }

    if (!first_iteration_has_been_done_) {
        auto mb_unitary_transformation_matrix = eigendecompositor_->BuildSubspectra(
            operators_to_calculate,
            derivatives_operators_to_calculate,
            number_of_block,
            subspace);

        eigenvectors_.push_back(mb_unitary_transformation_matrix);
        assert(eigenvectors_.size() == number_of_block + 1);

        const auto& energy_subspectrum =
            eigendecompositor_->getSpectrum(common::Energy)->get().blocks.at(number_of_block);

        auto raw_spectrum = energy_subspectrum.raw_data->multiply_by(1);
        current_energy_spectrum_.blocks.emplace_back(
            std::move(raw_spectrum),
            energy_subspectrum.properties);

        assert(current_energy_spectrum_.blocks.size() == number_of_block + 1);

        if (current_energy_derivative_spectrum_.has_value()) {
            auto raw_subspectrum_derivative =
                energy_subspectrum.raw_data->multiply_by(1 / initial_value_of_symbol_);
            current_energy_derivative_spectrum_.value().blocks.emplace_back(
                std::move(raw_subspectrum_derivative),
                energy_subspectrum.properties);
            assert(
                current_energy_derivative_spectrum_.value().blocks.size() == number_of_block + 1);
        }
    } else {
        double current_value_of_symbol = currentValueGetter_();
        double multiplier = current_value_of_symbol / initial_value_of_symbol_;
        const auto& energy_subspectrum =
            eigendecompositor_->getSpectrum(common::Energy)->get().blocks.at(number_of_block);
        auto raw_subspectrum_energy = energy_subspectrum.raw_data->multiply_by(multiplier);
        current_energy_spectrum_.blocks.emplace_back(
            std::move(raw_subspectrum_energy),
            energy_subspectrum.properties);
        assert(current_energy_spectrum_.blocks.size() == number_of_block + 1);

        // delete energy operator, because we just have implicitly calculated energy
        std::erase_if(operators_to_calculate, [](auto p) { return p.first == common::Energy; });
    }
    return eigenvectors_.at(number_of_block);
}

std::optional<std::reference_wrapper<const Spectrum>>
OneSymbolInHamiltonianEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::Energy) {
        return current_energy_spectrum_;
    }
    return eigendecompositor_->getSpectrum(quantity_enum);
}

std::optional<std::reference_wrapper<const Matrix>>
OneSymbolInHamiltonianEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::Energy) {
        return std::nullopt;
    }
    return eigendecompositor_->getMatrix(quantity_enum);
}

std::optional<std::reference_wrapper<const Spectrum>>
OneSymbolInHamiltonianEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (quantity_enum == common::Energy && current_energy_derivative_spectrum_.has_value()) {
        return current_energy_derivative_spectrum_.value();
    }
    return eigendecompositor_->getSpectrumDerivative(quantity_enum, symbol_name);
}

std::optional<std::reference_wrapper<const Matrix>>
OneSymbolInHamiltonianEigendecompositor::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (quantity_enum == common::Energy) {
        return std::nullopt;
    }
    return eigendecompositor_->getMatrixDerivative(quantity_enum, symbol_name);
}

void OneSymbolInHamiltonianEigendecompositor::initialize() {
    current_energy_spectrum_.blocks.clear();
    if (!first_iteration_has_been_done_) {
        eigendecompositor_->initialize();
    }
}

void OneSymbolInHamiltonianEigendecompositor::finalize() {
    if (!first_iteration_has_been_done_) {
        eigendecompositor_->finalize();
    }
    first_iteration_has_been_done_ = true;
}

}  // namespace eigendecompositor