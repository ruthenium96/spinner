#include "ImplicitSSquareEigendecompositor.h"

#include <cassert>
#include <utility>

namespace eigendecompositor {

ImplicitSSquareEigendecompositor::ImplicitSSquareEigendecompositor(
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
    quantum::linear_algebra::FactoriesList factories_list) :
    eigendecompositor_(std::move(eigendecompositor)),
    factories_list_(std::move(factories_list)) {
    s_square_implicit_spectrum_ = Spectrum();
}

std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
ImplicitSSquareEigendecompositor::BuildSubspectra(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>& operators_,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_,
    size_t number_of_block,
    const space::Subspace& subspace) {
    if (operators_.contains(common::S_total_squared)) {
        throw std::invalid_argument(
            "Explicit S2 operator passed to ImplicitSSquareEigendecompositor");
    }

    auto mb_unitary_transformation_matrix = eigendecompositor_->BuildSubspectra(
        operators_,
        derivatives_operators_,
        number_of_block,
        subspace);

    if (!first_iteration_has_been_done_) {
        // todo: we use energy_subspectrum only because of size_of_subspace
        //  can we implement calculation of size_of_subspace in Subspace class?
        const auto& energy_subspectrum =
            eigendecompositor_->getSpectrum(common::Energy)->get().blocks.at(number_of_block);

        double mult = subspace.properties.total_mult.value();
        double spin = (mult - 1) / 2.0;
        double s_squared_value = spin * (spin + 1);
        size_t size_of_subspace = energy_subspectrum.raw_data->size();

        auto raw_data = factories_list_.createVector();
        raw_data->add_identical_values(size_of_subspace, s_squared_value);
        s_square_implicit_spectrum_.blocks.emplace_back(std::move(raw_data), subspace.properties);
        assert(s_square_implicit_spectrum_.blocks.size() == number_of_block + 1);
    }

    return mb_unitary_transformation_matrix;
}

std::optional<std::reference_wrapper<const Spectrum>>
ImplicitSSquareEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::S_total_squared) {
        return s_square_implicit_spectrum_;
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

void ImplicitSSquareEigendecompositor::initialize() {
    eigendecompositor_->initialize();
}

void ImplicitSSquareEigendecompositor::finalize() {
    eigendecompositor_->finalize();
    first_iteration_has_been_done_ = true;
}

}  // namespace eigendecompositor