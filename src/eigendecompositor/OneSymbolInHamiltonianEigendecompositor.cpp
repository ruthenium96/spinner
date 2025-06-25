#include "OneSymbolInHamiltonianEigendecompositor.h"

#include <functional>
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

std::optional<OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
OneSymbolInHamiltonianEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {
    if (!first_iteration_has_been_done_) {
        auto mb_unitary_transformation_matrix = eigendecompositor_->BuildSubspectra(
            number_of_block,
            subspace);

        eigenvectors_[number_of_block] = mb_unitary_transformation_matrix;

        const auto energy_subspectrum = eigendecompositor_->getSubspectrum(common::Energy, number_of_block).value();

        current_energy_spectrum_[number_of_block] = transform_one_or_many(
            std::function([number_of_block](std::reference_wrapper<const Subspectrum> energy_subspectrum) {
                auto raw_spectrum = energy_subspectrum.get().raw_data->multiply_by(1);
                return Subspectrum(std::move(raw_spectrum), energy_subspectrum.get().properties);    
        }), energy_subspectrum);

        if (current_energy_derivative_spectrum_.has_value()) {
            current_energy_derivative_spectrum_.value()[number_of_block] = transform_one_or_many(
                std::function([number_of_block, this](std::reference_wrapper<const Subspectrum> energy_subspectrum) {
                    auto raw_subspectrum_derivative =
                        energy_subspectrum.get().raw_data->multiply_by(1 / initial_value_of_symbol_);
                    return Subspectrum(std::move(raw_subspectrum_derivative), energy_subspectrum.get().properties);
            }), energy_subspectrum);
        }
#ifndef NDEBUG
        const auto energy_submatrix = eigendecompositor_->getSubmatrix(common::Energy, number_of_block).value();
        current_energy_matrix_[number_of_block] = transform_one_or_many(
            std::function([number_of_block](std::reference_wrapper<const Submatrix> energy_submatrix) {
                auto raw_matrix = energy_submatrix.get().raw_data->multiply_by(1);
                return Submatrix(std::move(raw_matrix), energy_submatrix.get().properties);    
        }), energy_submatrix);
#endif
    } else {
        double current_value_of_symbol = currentValueGetter_();
        double multiplier = current_value_of_symbol / initial_value_of_symbol_;
        const auto energy_subspectrum = eigendecompositor_->getSubspectrum(common::Energy, number_of_block).value();
        current_energy_spectrum_[number_of_block] = transform_one_or_many(
            std::function([number_of_block, multiplier](std::reference_wrapper<const Subspectrum> energy_subspectrum) {
                auto raw_subspectrum = energy_subspectrum.get().raw_data->multiply_by(multiplier);
                return Subspectrum(std::move(raw_subspectrum), energy_subspectrum.get().properties);
        }), energy_subspectrum);
#ifndef NDEBUG
        const auto energy_submatrix = eigendecompositor_->getSubmatrix(common::Energy, number_of_block).value();
        current_energy_matrix_[number_of_block] = transform_one_or_many(
            std::function([number_of_block, multiplier](std::reference_wrapper<const Submatrix> energy_submatrix) {
                auto raw_submatrix = energy_submatrix.get().raw_data->multiply_by(multiplier);
                return Submatrix(std::move(raw_submatrix), energy_submatrix.get().properties);
        }), energy_submatrix);
#endif
    }
    return eigenvectors_.at(number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
OneSymbolInHamiltonianEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == common::Energy) {
        return copyRef<Subspectrum, std::reference_wrapper<const Subspectrum>>(current_energy_spectrum_[number_of_block]);
    }
    return eigendecompositor_->getSubspectrum(quantity_enum, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
OneSymbolInHamiltonianEigendecompositor::getSubmatrix(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == common::Energy) {
#ifndef NDEBUG
        return copyRef<Submatrix, std::reference_wrapper<const Submatrix>>(current_energy_matrix_[number_of_block]);
#else
        return std::nullopt;
#endif
    }
    return eigendecompositor_->getSubmatrix(quantity_enum, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
OneSymbolInHamiltonianEigendecompositor::getSubspectrumDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name, size_t number_of_block) const {
    if (quantity_enum == common::Energy && current_energy_derivative_spectrum_.has_value()) {
        return copyRef<Subspectrum, std::reference_wrapper<const Subspectrum>>(current_energy_derivative_spectrum_.value()[number_of_block]);
    }
    return eigendecompositor_->getSubspectrumDerivative(quantity_enum, symbol_name, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
OneSymbolInHamiltonianEigendecompositor::getSubmatrixDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name, size_t number_of_block) const {
    if (quantity_enum == common::Energy) {
        return std::nullopt;
    }
    return eigendecompositor_->getSubmatrixDerivative(quantity_enum, symbol_name, number_of_block);
}

void OneSymbolInHamiltonianEigendecompositor::initialize(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    uint32_t number_of_subspaces) {
    current_energy_spectrum_.clear();
    current_energy_spectrum_.resize(number_of_subspaces);

    size_t number_of_all_derivatives = derivatives_operators_to_calculate.size();
    // delete energy derivative operator, because we are about to calculate corresponding spectrum
    std::erase_if(derivatives_operators_to_calculate, [](const auto& p) {
        return p.first.first == common::Energy;
    });
    size_t number_of_non_energy_derivatives = derivatives_operators_to_calculate.size();
    size_t number_of_energy_derivatives =
        number_of_all_derivatives - number_of_non_energy_derivatives;

    if (number_of_energy_derivatives >= 2) {
        throw std::invalid_argument(
            "More than two energy derivatives passed to OneParameterEigendecompositor");
    }

    if (!first_iteration_has_been_done_) {
        eigendecompositor_->initialize(
            operators_to_calculate,
            derivatives_operators_to_calculate,
            number_of_subspaces);
        current_energy_derivative_spectrum_ = std::vector<OneOrMany<Subspectrum>>();
        current_energy_derivative_spectrum_->resize(number_of_subspaces);

        eigenvectors_.resize(number_of_subspaces);
    } else {
        // delete energy operator, because we are going to implicitly calculate energy
        std::erase_if(operators_to_calculate, [](const auto& p) { return p.first == common::Energy; });
    }
}

void OneSymbolInHamiltonianEigendecompositor::finalize() {
    if (!first_iteration_has_been_done_) {
        eigendecompositor_->finalize();
    }
    first_iteration_has_been_done_ = true;
}

}  // namespace eigendecompositor