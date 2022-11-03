#include "Model.h"

#include <cassert>
#include <utility>

#include "src/model/operators/terms/ScalarProductTerm.h"
#include "src/model/operators/terms/SzSzOneCenterTerm.h"
#include "src/model/operators/terms/SzSzTwoCenterTerm.h"

namespace model {
Model::Model(ModelInput modelInput) :
    modelInput_(std::move(modelInput)),
    converter_(modelInput_.getMults()) {
    energy_operator = operators::Operator();
    // TODO: this strange check need only because some tests do not initialize g factors,
    //  but want to calculate S^2 values. Fix it.
    if (getSymbols().isGFactorInitialized() && !getSymbols().isAllGFactorsEqual()) {
        InitializeGSzSquared();
        // TODO: when there is no Sz <-> -Sz symmetry, also \sum g_aS_{az} required
    } else {
        InitializeSSquared();
    }

    if (getSymbols().isZFSInitialized()) {
        InitializeZeroFieldSplitting();
    }

    if (getSymbols().isIsotropicExchangeInitialized()) {
        InitializeIsotropicExchange();
    }
}

void Model::InitializeDerivatives() {
    if (!getSymbols().isAllGFactorsEqual()) {
        InitializeGSzSquaredDerivatives();
    }

    if (getSymbols().isIsotropicExchangeInitialized()) {
        InitializeIsotropicExchangeDerivatives();
    }
}

void Model::InitializeSSquared() {
    if (operators_history_.s_squared) {
        return;
    }

    s_squared_operator = operators::Operator::s_squared(converter_);

    operators_history_.s_squared = true;
}

void Model::InitializeGSzSquared() {
    if (operators_history_.g_sz_squared) {
        return;
    }

    g_sz_squared_operator = operators::Operator::g_sz_squared(
        converter_,
        getSymbols().getGGParameters().first,
        getSymbols().getGGParameters().second);

    operators_history_.g_sz_squared = true;
}

void Model::InitializeIsotropicExchange() {
    if (operators_history_.isotropic_exchange_in_hamiltonian) {
        return;
    }
    energy_operator.getTwoCenterTerms().emplace_back(
        std::make_unique<const operators::ScalarProductTerm>(
            converter_,
            getSymbols().getIsotropicExchangeParameters()));
    operators_history_.isotropic_exchange_in_hamiltonian = true;
}

void Model::InitializeZeroFieldSplitting() {
    if (operators_history_.zfs_in_hamiltonian) {
        return;
    }
    energy_operator.getOneCenterTerms().emplace_back(
        std::make_unique<const operators::SzSzOneCenterTerm>(
            converter_,
            getSymbols().getZFSParameters().first));
    operators_history_.zfs_in_hamiltonian = true;
}

void Model::InitializeIsotropicExchangeDerivatives() {
    if (operators_history_.isotropic_exchange_derivatives) {
        return;
    }

    for (const auto& symbol : getSymbols().getChangeableNames(symbols::SymbolTypeEnum::J)) {
        operators::Operator operator_derivative = operators::Operator();
        operator_derivative.getTwoCenterTerms().emplace_back(
            std::make_unique<const operators::ScalarProductTerm>(
                converter_,
                getSymbols().constructIsotropicExchangeDerivativeParameters(symbol)));
        derivatives_map_[{common::Energy, symbol}] = std::move(operator_derivative);
    }

    operators_history_.isotropic_exchange_derivatives = true;
}

void Model::InitializeGSzSquaredDerivatives() {
    if (operators_history_.g_sz_squared_derivatives) {
        return;
    }

    for (const auto& symbol : getSymbols().getChangeableNames(symbols::SymbolTypeEnum::g_factor)) {
        operators::Operator operator_derivative = operators::Operator();
        auto pair_of_parameters = getSymbols().constructGGDerivativeParameters(symbol);
        operator_derivative.getOneCenterTerms().emplace_back(
            std::make_unique<operators::SzSzOneCenterTerm>(converter_, pair_of_parameters.first));
        operator_derivative.getTwoCenterTerms().emplace_back(
            std::make_unique<const operators::SzSzTwoCenterTerm>(
                converter_,
                pair_of_parameters.second));
        derivatives_map_[{common::gSz_total_squared, symbol}] = std::move(operator_derivative);
    }

    operators_history_.g_sz_squared_derivatives = true;
}

const lexicographic::IndexConverter& Model::getIndexConverter() const {
    return converter_;
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

bool Model::is_zero_field_splitting_initialized() const {
    return operators_history_.zfs_in_hamiltonian;
}

const symbols::Symbols& Model::getSymbols() const {
    return modelInput_.getSymbols();
}

symbols::Symbols& Model::getSymbols() {
    return modelInput_.getSymbols();
}
}  // namespace model