#include "ImplicitQuantityEigendecompositor.h"

#include <sstream>
#include <stdexcept>
#include <utility>
#include "magic_enum.hpp"
#include "src/common/Quantity.h"

namespace eigendecompositor {

ImplicitQuantityEigendecompositor::ImplicitQuantityEigendecompositor(
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
    quantum::linear_algebra::FactoriesList factories_list,
    common::QuantityEnum quantity_implicit_enum,
    uint32_t max_ntz_proj) :
    eigendecompositor_(std::move(eigendecompositor)),
    max_ntz_proj_(max_ntz_proj),
    factories_list_(std::move(factories_list)),
    quantity_implicit_enum_(quantity_implicit_enum) {
    if (quantity_implicit_enum_ != common::S_total_squared 
        && quantity_implicit_enum_ != common::M_total_squared) {
        std::ostringstream os; 
        os << "ImplicitQuantityEigendecompositor cannot construct ";
        os << magic_enum::enum_name(quantity_implicit_enum_);
        os << " quantity\n";
        throw std::invalid_argument(os.str());
    }
    quantity_implicit_spectrum_ = Spectrum();
}

std::optional<OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
ImplicitQuantityEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {

    auto mb_unitary_transformation_matrix = eigendecompositor_->BuildSubspectra(
        number_of_block,
        subspace);

    if (!first_iteration_has_been_done_) {
        double value = calculate_value(subspace.properties);

        size_t size_of_subspectrum = getSubspectrumSize(common::Energy, number_of_block);

        auto raw_data = factories_list_.createVector();
        raw_data->add_identical_values(size_of_subspectrum, value);
        quantity_implicit_spectrum_.blocks[number_of_block] =
            Subspectrum(std::move(raw_data), subspace.properties);
    }

    return mb_unitary_transformation_matrix;
}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
ImplicitQuantityEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == quantity_implicit_enum_) {
        return quantity_implicit_spectrum_.blocks[number_of_block];
    }
    return eigendecompositor_->getSubspectrum(quantity_enum, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
ImplicitQuantityEigendecompositor::getSubmatrix(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == quantity_implicit_enum_) {
        return std::nullopt;
    }
    return eigendecompositor_->getSubmatrix(quantity_enum, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
ImplicitQuantityEigendecompositor::getSubspectrumDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name, size_t number_of_block) const {
    if (quantity_enum == quantity_implicit_enum_) {
        return std::nullopt;
    }
    return eigendecompositor_->getSubspectrumDerivative(quantity_enum, symbol_name, number_of_block);
}

std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
ImplicitQuantityEigendecompositor::getSubmatrixDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name, size_t number_of_block) const {
    if (quantity_enum == quantity_implicit_enum_) {
        return std::nullopt;
    }
    return eigendecompositor_->getSubmatrixDerivative(quantity_enum, symbol_name, number_of_block);
}

OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>
ImplicitQuantityEigendecompositor::getWeightsOfBlockStates(size_t number_of_block) const {
    return eigendecompositor_->getWeightsOfBlockStates(number_of_block);
}

void ImplicitQuantityEigendecompositor::initialize(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    uint32_t number_of_subspaces) {
    if (operators_to_calculate.contains(common::S_total_squared)) {
        throw std::invalid_argument(
            "Explicit S^2 operator passed to ImplicitSSquareEigendecompositor");
    }
    if (operators_to_calculate.contains(common::M_total_squared)) {
        throw std::invalid_argument(
            "Explicit M^2 operator passed to ImplicitSSquareEigendecompositor");
    }
    if (!first_iteration_has_been_done_) {
        quantity_implicit_spectrum_.blocks.resize(number_of_subspaces);
    }
    eigendecompositor_->initialize(
        operators_to_calculate,
        derivatives_operators_to_calculate,
        number_of_subspaces);
}

void ImplicitQuantityEigendecompositor::finalize() {
    eigendecompositor_->finalize();
    first_iteration_has_been_done_ = true;
}

double ImplicitQuantityEigendecompositor::calculate_value(const BlockProperties& properties) const {
    if (quantity_implicit_enum_ == common::S_total_squared) {
        double mult = properties.total_mult.value();
        double spin = (mult - 1) / 2.0;
        return spin * (spin + 1);    
    }
    if (quantity_implicit_enum_ == common::M_total_squared) {
        auto max_spin = ((double)max_ntz_proj_ - 1.0) / 2.0;
        double total_projection = properties.n_proj.value() - max_spin;
        return total_projection * total_projection;    
    }
}


}  // namespace eigendecompositor