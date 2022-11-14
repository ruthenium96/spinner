#include "Model.h"

#include <cassert>
#include <utility>

#include "src/model/operators/terms/ScalarProductTerm.h"
#include "src/model/operators/terms/SzSzOneCenterTerm.h"
#include "src/model/operators/terms/SzSzTwoCenterTerm.h"

namespace model {
Model::Model(ModelInput modelInput) :
    numericalWorker_(modelInput.modifySymbolicWorker(), modelInput.getMults().size()),
    converter_(modelInput.getMults()) {
    energy_operator = operators::Operator();
    // TODO: this strange check need only because some tests do not initialize g factors,
    //  but want to calculate S^2 values. Fix it.
    if (getSymbolicWorker().isGFactorInitialized()
        && (!getSymbolicWorker().isAllGFactorsEqual() || getSymbolicWorker().isZFSInitialized())) {
        InitializeGSzSquared();
        // TODO: when there is no Sz <-> -Sz symmetry, also \sum g_aS_{az} required
    } else {
        InitializeSSquared();
    }

    if (getSymbolicWorker().isZFSInitialized()) {
        InitializeZeroFieldSplitting();
    }

    if (getSymbolicWorker().isIsotropicExchangeInitialized()) {
        InitializeIsotropicExchange();
    }
}

void Model::InitializeDerivatives() {
    if (!getSymbolicWorker().isAllGFactorsEqual()) {
        InitializeGSzSquaredDerivatives();
    }

    if (getSymbolicWorker().isIsotropicExchangeInitialized()) {
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
        getNumericalWorker().getGGParameters().first,
        getNumericalWorker().getGGParameters().second);

    operators_history_.g_sz_squared = true;
}

void Model::InitializeIsotropicExchange() {
    if (operators_history_.isotropic_exchange_in_hamiltonian) {
        return;
    }
    energy_operator.getTwoCenterTerms().emplace_back(
        std::make_unique<const operators::ScalarProductTerm>(
            converter_,
            getNumericalWorker().getIsotropicExchangeParameters()));
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

    for (const auto& symbol : getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::J)) {
        operators::Operator operator_derivative = operators::Operator();
        operator_derivative.getTwoCenterTerms().emplace_back(
            std::make_unique<const operators::ScalarProductTerm>(
                converter_,
                getNumericalWorker().constructIsotropicExchangeDerivativeParameters(symbol)));
        derivatives_map_[{common::Energy, symbol}] = std::move(operator_derivative);
    }

    operators_history_.isotropic_exchange_derivatives = true;
}

void Model::InitializeGSzSquaredDerivatives() {
    if (operators_history_.g_sz_squared_derivatives) {
        return;
    }

    for (const auto& symbol :
         getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::g_factor)) {
        operators::Operator operator_derivative = operators::Operator();
        auto pair_of_parameters = getNumericalWorker().constructGGDerivativeParameters(symbol);
        operator_derivative.getOneCenterTerms().emplace_back(
            std::make_unique<operators::SzSzOneCenterTerm>(converter_, pair_of_parameters.first));
        operator_derivative.getTwoCenterTerms().emplace_back(
            std::make_unique<const operators::SzSzTwoCenterTerm>(
                converter_,
                pair_of_parameters.second,
                2));
        // this two from summation in Submatrix: \sum_{a=1}^N \sum_{b=a+1}^N
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

const symbols::SymbolicWorker& Model::getSymbolicWorker() const {
    return numericalWorker_.getSymbolicWorker();
}

symbols::NumericalWorker& Model::getNumericalWorker() {
    return numericalWorker_;
}

const symbols::NumericalWorker& Model::getNumericalWorker() const {
    return numericalWorker_;
}
}  // namespace model