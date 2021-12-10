#include "Symbols.h"

#include <cmath>
#include <utility>

namespace symbols {
Symbols::Symbols(size_t number_of_spins) : number_of_spins_(number_of_spins) {}

void Symbols::assignSymbolToIsotropicExchange(
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

void Symbols::assignSymbolToGFactor(const SymbolName& symbol_name, size_t center_a) {
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
}

std::shared_ptr<const DenseVector> Symbols::getGFactorParameters() const {
    if (numeric_g_factors_ == nullptr) {
        throw std::invalid_argument("g-factors were not initialized");
    }
    return numeric_g_factors_;
}

bool Symbols::symmetry_consistence(const group::Group& group) const {
    // isotropic exchange parameters:
    if (!symbolic_isotropic_exchanges_.empty()) {
        for (const auto& element : group.elements_) {
            std::vector<std::vector<SymbolName>> permutated_symbols(
                element.size(),
                std::vector<SymbolName>(element.size()));
            for (size_t i = 0; i < element.size(); ++i) {
                for (size_t j = 0; j < element.size(); ++j) {
                    permutated_symbols[element[i]][element[j]] =
                        symbolic_isotropic_exchanges_[i][j];
                }
            }
            if (permutated_symbols != symbolic_isotropic_exchanges_) {
                return false;
            }
        }
    }
    if (!symbolic_g_factors_.empty()) {
        for (const auto& element : group.elements_) {
            std::vector<SymbolName> permutated_symbols(element.size());
            for (size_t i = 0; i < element.size(); ++i) {
                permutated_symbols[element[i]] = symbolic_g_factors_[i];
            }
            if (permutated_symbols != symbolic_g_factors_) {
                return false;
            }
        }
    }
    return true;
}

std::shared_ptr<const DenseMatrix>
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

    auto ptr_to_derivative = std::make_shared<DenseMatrix>();
    ptr_to_derivative->resize_with_nans(number_of_spins_, number_of_spins_);

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            if (symbolic_isotropic_exchanges_[i][j] != symbol_name) {
                value = NAN;
            } else {
                value = 1;
            }
            ptr_to_derivative->assign_to_position(value, i, j);
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
    if (symbolic_g_factors_.empty()) {
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

void Symbols::setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value) {
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
        // and other things
    }
}

void Symbols::updateIsotropicExchangeParameters() {
    if (numeric_isotropic_exchanges_ == nullptr) {
        numeric_isotropic_exchanges_ = std::make_shared<DenseMatrix>();
        numeric_isotropic_exchanges_->resize_with_nans(number_of_spins_, number_of_spins_);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            if (symbolic_isotropic_exchanges_[i][j].get_name().empty()) {
                value = NAN;
            } else {
                value = symbolsMap[symbolic_isotropic_exchanges_[i][j]].value;
            }
            numeric_isotropic_exchanges_->assign_to_position(value, i, j);
        }
    }
}

void Symbols::updateGFactorParameters() {
    if (numeric_g_factors_ == nullptr) {
        numeric_g_factors_ = std::make_shared<DenseVector>();
        numeric_g_factors_->resize(number_of_spins_);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        double value = symbolsMap[symbolic_g_factors_[i]].value;
        numeric_g_factors_->assign_to_position(value, i);
    }
}

std::shared_ptr<const DenseMatrix> Symbols::getIsotropicExchangeParameters() const {
    if (numeric_isotropic_exchanges_ == nullptr) {
        throw std::invalid_argument("Isotropic exchange interaction was not initialized");
    }
    return numeric_isotropic_exchanges_;
}

}  // namespace symbols