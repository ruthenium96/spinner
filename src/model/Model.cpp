#include "Model.h"

#include <memory>
#include <utility>

#include "src/common/Quantity.h"
#include "src/model/operators/Operator.h"

namespace model {
Model::Model(ModelInput modelInput, std::unique_ptr<operators::AbstractOperatorConstructor>&& operator_constructor) :
    numericalWorker_(modelInput.getSymbolicWorker(), modelInput.getMults().size()),
    operator_constructor_(std::move(operator_constructor)) {
    operators_map_[common::Energy] = std::make_shared<operators::Operator>();
    operators_map_[common::S_total_squared] = operator_constructor_->constructSSquared();
    
    if (getSymbolicWorker().isGFactorInitialized()) {
        operators_map_[common::gSz_total_squared] = operator_constructor_->constructGSzSquaredLike(
            getNumericalWorker().getGGParameters().first,
            getNumericalWorker().getGGParameters().second);

        for (const auto& symbol :
            getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::g_factor)) {
            auto pair_of_parameters = getNumericalWorker().constructGGDerivativeParameters(symbol);
            derivatives_map_[{common::gSz_total_squared, symbol}] = operator_constructor_->constructGSzSquaredLike(
                pair_of_parameters.first,
                pair_of_parameters.second
            );
        }
    }

    if (getSymbolicWorker().isZFSInitialized()) {
        operator_constructor_->emplaceZeroFieldSplittingLike(
            operators_map_[common::Energy], 
            getNumericalWorker().getZFSParameters().first);

        for (const auto& symbol : getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::D)) {
            auto operator_derivative = std::make_shared<operators::Operator>();
            auto derivative_parameters = getNumericalWorker().constructZFSDerivativeParameters(symbol);
            operator_constructor_->emplaceZeroFieldSplittingLike(operator_derivative, derivative_parameters);
            derivatives_map_[{common::Energy, symbol}] = operator_derivative;
        }
    }

    if (getSymbolicWorker().isIsotropicExchangeInitialized()) {
        operator_constructor_->emplaceIsotropicExchangeLike(
            operators_map_[common::Energy],
            getNumericalWorker().getIsotropicExchangeParameters());

        for (const auto& symbol : getSymbolicWorker().getChangeableNames(symbols::SymbolTypeEnum::J)) {
            auto operator_derivative = std::make_shared<operators::Operator>();
            operator_constructor_->emplaceIsotropicExchangeLike(
                operator_derivative, 
                getNumericalWorker().constructIsotropicExchangeDerivativeParameters(symbol));
            derivatives_map_[{common::Energy, symbol}] = operator_derivative;
        }
    }
}

std::optional<std::shared_ptr<const operators::Operator>>
Model::getOperator(common::QuantityEnum quantity_enum) const {
    if (operators_map_.contains(quantity_enum)) {
        return operators_map_.at(quantity_enum);
    }
    return std::nullopt;
}

const std::map<common::QuantityEnum, std::shared_ptr<operators::Operator>>& Model::getOperators() {
    return operators_map_;
}

std::optional<std::shared_ptr<const operators::Operator>> Model::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    const symbols::SymbolName& symbol) const {
    if (derivatives_map_.contains({quantity_enum, symbol})) {
        return derivatives_map_.at({quantity_enum, symbol});
    }
    return std::nullopt;
}

const std::
    map<std::pair<common::QuantityEnum, symbols::SymbolName>, std::shared_ptr<operators::Operator>>&
    Model::getOperatorDerivatives() {
    return derivatives_map_;
}


bool Model::is_g_sz_squared_initialized() const {
    return getSymbolicWorker().isGFactorInitialized();
}

bool Model::is_g_sz_squared_derivatives_initialized() const {
    return getSymbolicWorker().isGFactorInitialized();
}

bool Model::gFactorsAreAllNoneOrAreTheSame() const {
    return !getSymbolicWorker().isGFactorInitialized() ||
    (getSymbolicWorker().isGFactorInitialized() && getSymbolicWorker().isAllGFactorsEqual());
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

void Model::setNewValueToChangeableSymbol(
    const symbols::SymbolName& symbol_name,
    double new_value) {
    numericalWorker_.setNewValueToChangeableSymbol(symbol_name, new_value);
}
}  // namespace model