#include "Symbols.h"

#include <cmath>
#include <utility>

namespace model::symbols {
Symbols::Symbols(size_t number_of_spins) : number_of_spins_(number_of_spins) {}

Symbols& Symbols::assignSymbolToIsotropicExchange(
    const SymbolName& symbol_name,
    size_t center_a,
    size_t center_b) {
    if (symbolsMap.find(symbol_name) == symbolsMap.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has not been initialized");
    }
    if (symbolic_isotropic_exchanges_.empty()) {
        symbolic_isotropic_exchanges_.resize(
            number_of_spins_,
            std::vector<SymbolName>(number_of_spins_));
    }
    if (center_b == center_a) {
        throw std::invalid_argument("Isotropic exchange takes place between different centers");
    }
    if (!symbolic_isotropic_exchanges_[center_a][center_b].get_name().empty()) {
        throw std::invalid_argument("This parameter has been already specified");
    }

    if (symbolsMap[symbol_name].type_enum != SymbolTypeEnum::not_specified
        && symbolsMap[symbol_name].type_enum != SymbolTypeEnum::J) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been already specified as not J parameter");
    }
    if (symbolsMap[symbol_name].type_enum == SymbolTypeEnum::not_specified) {
        symbolsMap[symbol_name].type_enum = SymbolTypeEnum::J;
    }

    symbolic_isotropic_exchanges_[center_a][center_b] = symbol_name;
    symbolic_isotropic_exchanges_[center_b][center_a] = symbol_name;

    updateIsotropicExchangeParameters();

    return *this;
}

Symbols& Symbols::assignSymbolToGFactor(const SymbolName& symbol_name, size_t center_a) {
    if (symbolsMap.find(symbol_name) == symbolsMap.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has not been initialized");
    }

    if (symbolsMap[symbol_name].type_enum != SymbolTypeEnum::not_specified
        && symbolsMap[symbol_name].type_enum != SymbolTypeEnum::g_factor) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been already specified as not g factor parameter");
    }
    if (symbolsMap[symbol_name].type_enum == SymbolTypeEnum::not_specified) {
        symbolsMap[symbol_name].type_enum = SymbolTypeEnum::g_factor;
    }

    if (symbolic_g_factors_.empty()) {
        symbolic_g_factors_.resize(number_of_spins_);
    }

    symbolic_g_factors_[center_a] = symbol_name;

    updateGFactorParameters();

    return *this;
}

Symbols& Symbols::assignSymbolToTheta(const SymbolName& symbol_name) {
    if (symbolsMap.find(symbol_name) == symbolsMap.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has not been initialized");
    }

    if (symbolsMap[symbol_name].type_enum != SymbolTypeEnum::not_specified
        && symbolsMap[symbol_name].type_enum != SymbolTypeEnum::Theta) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been already specified as not Theta parameter");
    }
    if (symbolsMap[symbol_name].type_enum == SymbolTypeEnum::not_specified) {
        symbolsMap[symbol_name].type_enum = SymbolTypeEnum::Theta;
    }

    symbolic_Theta_ = symbol_name;

    updateThetaParameter();

    return *this;
}

SymbolName Symbols::addSymbol(
    const std::string& name_string,
    double initial_value,
    bool is_changeable,
    SymbolTypeEnum type_enum) {
    SymbolName symbol_name(name_string);
    if (symbolsMap.find(symbol_name) != symbolsMap.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has been already initialized");
    }
    symbolsMap[symbol_name] = {initial_value, is_changeable, type_enum};
    return symbol_name;
}

SymbolName
Symbols::addSymbol(const std::string& name_string, double initial_value, bool is_changeable) {
    return addSymbol(name_string, initial_value, is_changeable, SymbolTypeEnum::not_specified);
}

SymbolName Symbols::addSymbol(const std::string& name_string, double initial_value) {
    return addSymbol(name_string, initial_value, true);
}

std::shared_ptr<const OneDNumericalParameters<double>> Symbols::getGFactorParameters() const {
    if (numeric_g_factors_ == nullptr) {
        throw std::invalid_argument("g-factors were not initialized");
    }
    return numeric_g_factors_;
}

std::shared_ptr<const TwoDNumericalParameters<double>>
Symbols::getGGFactorProductParameters() const {
    if (numeric_g_g_factors_product_ == nullptr) {
        throw std::invalid_argument("g-factors were not initialized");
    }
    return numeric_g_g_factors_product_;
}

std::shared_ptr<const TwoDNumericalParameters<double>>
Symbols::constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name) {
    if (symbolsMap.find(symbol_name) == symbolsMap.end()) {
        throw std::invalid_argument(symbol_name.get_name() + " name has not been initialized");
    }

    if (symbolsMap[symbol_name].type_enum != SymbolTypeEnum::J) {
        throw std::invalid_argument(
            symbol_name.get_name() + " has been specified as not J parameter");
    }

    if (symbolic_isotropic_exchanges_.empty()) {
        throw std::length_error("Isotropic exchange parameters has not been initialized");
    }

    auto ptr_to_derivative =
        std::make_shared<TwoDNumericalParameters<double>>(number_of_spins_, NAN);

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            if (symbolic_isotropic_exchanges_[i][j] != symbol_name) {
                value = NAN;
            } else {
                value = 1;
            }
            ptr_to_derivative->at(i, j) = value;
        }
    }

    numeric_isotropic_exchange_derivatives_.push_back(ptr_to_derivative);

    return ptr_to_derivative;
}

std::vector<SymbolName> Symbols::getChangeableNames(SymbolTypeEnum type_enum) const {
    std::vector<SymbolName> answer;
    for (const auto& [symbol_name, symbol_data] : symbolsMap) {
        if (symbol_data.is_changeable && symbol_data.type_enum == type_enum) {
            answer.push_back(symbol_name);
        }
    }
    return answer;
}

bool Symbols::isAllGFactorsEqual() const {
    if (!isGFactorInitialized()) {
        throw std::length_error("g factor parameters have not been initialized");
    }
    for (size_t i = 1; i < symbolic_g_factors_.size(); ++i) {
        if (symbolic_g_factors_[i] != symbolic_g_factors_[0]) {
            return false;
        }
    }
    return true;
}

std::vector<SymbolName> Symbols::getChangeableNames() const {
    std::vector<SymbolName> answer;
    for (const auto& [symbol_name, symbol_data] : symbolsMap) {
        if (symbol_data.is_changeable) {
            answer.push_back(symbol_name);
        }
    }
    return answer;
}

double Symbols::getValueOfName(const SymbolName& symbol_name) const {
    return symbolsMap.at(symbol_name).value;
}

Symbols& Symbols::setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value) {
    SymbolData symbol_data = symbolsMap[symbol_name];
    if (!symbol_data.is_changeable) {
        throw std::invalid_argument("Cannot change value of unchangeable symbol");
    }
    symbol_data.value = new_value;
    symbolsMap[symbol_name] = symbol_data;

    if (symbol_data.type_enum == symbols::SymbolTypeEnum::J) {
        updateIsotropicExchangeParameters();
        // NB: numeric_isotropic_exchange_derivatives_ does not change when J changes.
    } else if (symbol_data.type_enum == symbols::SymbolTypeEnum::g_factor) {
        updateGFactorParameters();
    } else if (symbol_data.type_enum == symbols::SymbolTypeEnum::Theta) {
        updateThetaParameter();
    }
    // and other things

    return *this;
}

void Symbols::updateIsotropicExchangeParameters() {
    if (numeric_isotropic_exchanges_ == nullptr) {
        numeric_isotropic_exchanges_ =
            std::make_shared<TwoDNumericalParameters<double>>(number_of_spins_, NAN);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            if (symbolic_isotropic_exchanges_[i][j].get_name().empty()) {
                value = NAN;
            } else {
                value = symbolsMap[symbolic_isotropic_exchanges_[i][j]].value;
            }
            numeric_isotropic_exchanges_->at(i, j) = value;
        }
    }
}

void Symbols::updateGFactorParameters() {
    if (numeric_g_factors_ == nullptr) {
        numeric_g_factors_ =
            std::make_shared<OneDNumericalParameters<double>>(number_of_spins_, NAN);
    }
    if (numeric_g_g_factors_product_ == nullptr) {
        numeric_g_g_factors_product_ =
            std::make_shared<TwoDNumericalParameters<double>>(number_of_spins_, NAN);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        double value = symbolsMap[symbolic_g_factors_[i]].value;
        numeric_g_factors_->at(i) = value;
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value =
                symbolsMap[symbolic_g_factors_[i]].value * symbolsMap[symbolic_g_factors_[j]].value;
            numeric_g_g_factors_product_->at(i, j) = value;
        }
    }
}

void Symbols::updateThetaParameter() {
    if (numeric_Theta_ == nullptr) {
        numeric_Theta_ = std::make_shared<double>(NAN);
    }

    double value = symbolsMap[symbolic_Theta_.value()].value;
    *numeric_Theta_ = value;
}

std::shared_ptr<const TwoDNumericalParameters<double>>
Symbols::getIsotropicExchangeParameters() const {
    // TODO: refactor it:
    if (numeric_isotropic_exchanges_ == nullptr) {
        throw std::invalid_argument("Isotropic exchange interaction was not initialized");
    }
    return numeric_isotropic_exchanges_;
}

std::shared_ptr<const double> Symbols::getThetaParameter() const {
    if (numeric_Theta_ == nullptr) {
        throw std::invalid_argument("Theta was not initialized");
    }
    return numeric_Theta_;
}

bool Symbols::isIsotropicExchangeInitialized() const {
    return !symbolic_isotropic_exchanges_.empty();
}

bool Symbols::isGFactorInitialized() const {
    return !symbolic_g_factors_.empty();
}

bool Symbols::isThetaInitialized() const {
    return symbolic_Theta_.has_value();
}

SymbolName Symbols::getIsotropicExchangeSymbolName(size_t i, size_t j) const {
    return symbolic_isotropic_exchanges_[i][j];
}

SymbolName Symbols::getGFactorSymbolName(size_t i) const {
    return symbolic_g_factors_[i];
}

SymbolName Symbols::getThetaSymbolName() const {
    return symbolic_Theta_.value();
}

}  // namespace symbols