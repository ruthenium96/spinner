#include "Model.h"

#include <cassert>

#include "src/model/operators/terms/ScalarProductTerm.h"
#include "src/model/operators/terms/SzSzOneCenterTerm.h"
#include "src/model/operators/terms/SzSzTwoCenterTerm.h"

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

Model& Model::InitializeGSzSquared() {
    if (operators_history_.g_sz_squared) {
        return *this;
    }

    g_sz_squared_operator = operators::Operator::g_sz_squared(
        converter_,
        symbols_.getGGParameters().first,
        symbols_.getGGParameters().second);

    operators_history_.g_sz_squared = true;
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
        derivatives_map_[{common::Energy, symbol}] = std::move(operator_derivative);
    }

    operators_history_.isotropic_exchange_derivatives = true;

    return *this;
}

Model& Model::InitializeGSzSquaredDerivatives() {
    if (operators_history_.g_sz_squared_derivatives) {
        return *this;
    }

    for (const auto& symbol : symbols_.getChangeableNames(symbols::SymbolTypeEnum::g_factor)) {
        operators::Operator operator_derivative = operators::Operator();
        auto pair_of_parameters = symbols_.constructGGDerivativeParameters(symbol);
        operator_derivative.getOneCenterTerms().emplace_back(
            std::make_unique<operators::SzSzOneCenterTerm>(converter_, pair_of_parameters.first));
        operator_derivative.getTwoCenterTerms().emplace_back(
            std::make_unique<const operators::SzSzTwoCenterTerm>(
                converter_,
                pair_of_parameters.second));
        derivatives_map_[{common::gSz_total_squared, symbol}] = std::move(operator_derivative);
    }

    operators_history_.g_sz_squared_derivatives = true;

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
    } else if (quantity_enum == common::QuantityEnum::gSz_total_squared) {
        return g_sz_squared_operator.value();
    }
    assert(0);
}

const operators::Operator& Model::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    const symbols::SymbolName& symbol) const {
    return derivatives_map_.at({quantity_enum, symbol});
}

bool Model::is_s_squared_initialized() const {
    return operators_history_.s_squared;
}

bool Model::is_isotropic_exchange_derivatives_initialized() const {
    return operators_history_.isotropic_exchange_derivatives;
}

bool Model::is_g_sz_squared_initialized() const {
    return operators_history_.g_sz_squared;
}

bool Model::is_g_sz_squared_derivatives_initialized() const {
    return operators_history_.g_sz_squared_derivatives;
}

}  // namespace model