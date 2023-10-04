#include "ImplicitSSquareEigendecompositor.h"

#include <utility>

namespace eigendecompositor {

ImplicitSSquareEigendecompositor::ImplicitSSquareEigendecompositor(
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
    quantum::linear_algebra::FactoriesList factories_list) :
    eigendecompositor_(std::move(eigendecompositor)),
    factories_list_(std::move(factories_list)) {}

void ImplicitSSquareEigendecompositor::BuildSpectra(
    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_,
    const space::Space& space) {
    if (operators_.contains(common::S_total_squared)) {
        throw std::invalid_argument(
            "Explicit S2 operator passed to ImplicitSSquareEigendecompositor");
    }

    eigendecompositor_->BuildSpectra(operators_, derivatives_operators_, space);

    for (const auto& subspectrum_energy : getSpectrum(common::Energy)->get().blocks) {
        double mult = subspectrum_energy.properties.total_mult.value();
        double spin = (mult - 1) / 2.0;
        double s_squared_value = spin * (spin + 1);
        size_t size = subspectrum_energy.raw_data->size();

        auto raw_data = factories_list_.createVector();
        raw_data->add_identical_values(size, s_squared_value);
        s_square_implicit_.spectrum_.blocks.emplace_back(
            std::move(raw_data),
            subspectrum_energy.properties);
    }
}

std::optional<std::reference_wrapper<const Spectrum>>
ImplicitSSquareEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::S_total_squared) {
        return s_square_implicit_.spectrum_;
    }
    return eigendecompositor_->getSpectrum(quantity_enum);
}

std::optional<std::reference_wrapper<const Matrix>>
ImplicitSSquareEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::S_total_squared) {
        return std::nullopt;
    }
    return eigendecompositor_->getMatrix(quantity_enum);
}

std::optional<std::reference_wrapper<const Spectrum>>
ImplicitSSquareEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (quantity_enum == common::S_total_squared) {
        return std::nullopt;
    }
    return eigendecompositor_->getSpectrumDerivative(quantity_enum, symbol_name);
}

std::optional<std::reference_wrapper<const Matrix>>
ImplicitSSquareEigendecompositor::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol_name) const {
    if (quantity_enum == common::S_total_squared) {
        return std::nullopt;
    }
    return eigendecompositor_->getMatrixDerivative(quantity_enum, symbol_name);
}
}  // namespace eigendecompositor