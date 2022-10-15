#include "Model.h"

#include <src/model/operators/ScalarProductTerm.h>

#include <cassert>

namespace model {

Model::Model(const std::vector<int>& mults) : symbols_(mults.size()), converter_(mults) {
    energy_operator = operators::Operator();
}

Model& Model::InitializeSSquared() {
    if (operators_history_.s_squared) {
        return *this;
    }

    s_squared_operator = operators::Operator::s_squared(converter_);

    operators_history_.s_squared = true;
    return *this;
}

Model& Model::InitializeIsotropicExchange() {
    if (operators_history_.isotropic_exchange_in_hamiltonian) {
        return *this;
    }
    energy_operator.getTwoCenterTerms().emplace_back(
        std::make_unique<const operators::ScalarProductTerm>(
            converter_,
            symbols_.getIsotropicExchangeParameters()));
    operators_history_.isotropic_exchange_in_hamiltonian = true;
    return *this;
}

Model& Model::InitializeIsotropicExchangeDerivatives() {
    if (operators_history_.isotropic_exchange_derivatives) {
        return *this;
    }

    for (const auto& symbol : symbols_.getChangeableNames(symbols::SymbolTypeEnum::J)) {
        operators::Operator operator_derivative = operators::Operator();
        operator_derivative.getTwoCenterTerms().emplace_back(
            std::make_unique<const operators::ScalarProductTerm>(
                converter_,
                symbols_.constructIsotropicExchangeDerivativeParameters(symbol)));
        derivative_of_energy_wrt_exchange_parameters_operator[symbol] =
            std::move(operator_derivative);
    }

    operators_history_.isotropic_exchange_derivatives = true;

    return *this;
}

const lexicographic::IndexConverter& Model::getIndexConverter() const {
    return converter_;
}

const symbols::Symbols& Model::getSymbols() const {
    return symbols_;
}

symbols::Symbols& Model::getSymbols() {
    return symbols_;
}

const operators::Operator& Model::getOperator(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return energy_operator;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return s_squared_operator.value();
    }
    assert(0);
}

const operators::Operator& Model::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    symbols::SymbolTypeEnum symbol_type,
    const symbols::SymbolName& symbol) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        if (symbol_type == symbols::SymbolTypeEnum::J) {
            return derivative_of_energy_wrt_exchange_parameters_operator.at(symbol);
        }
    }
    assert(0);
}

bool Model::is_s_squared_initialized() const {
    return operators_history_.s_squared;
}

bool Model::is_isotropic_exchange_derivatives_initialized() const {
    return operators_history_.isotropic_exchange_derivatives;
}

}  // namespace model