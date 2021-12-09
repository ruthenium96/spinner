#include "Symbols.h"

#include <cmath>
#include <utility>

namespace symbols {
Symbols::Symbols(size_t number_of_spins) : number_of_spins_(number_of_spins) {}

std::shared_ptr<const DenseMatrix> symbols::Symbols::constructIsotropicExchangeParameters() {
    if (isotropic_exchange_parameters_names_.empty()) {
        throw std::length_error("Isotropic exchange parameters has not been initialized");
    }

    if (isotropic_exchange_parameters_values_ == nullptr) {
        isotropic_exchange_parameters_values_ = std::make_shared<DenseMatrix>();
        isotropic_exchange_parameters_values_->resize_with_nans(number_of_spins_, number_of_spins_);
    }

    updateIsotropicExchangeParameters();

    return isotropic_exchange_parameters_values_;
}

void Symbols::addIsotropicExchange(
    const SymbolName& symbol_name,
    size_t center_a,
    size_t center_b) {
    if (symbols_.find(symbol_name) == symbols_.end()) {
        throw std::invalid_argument(symbol_name.name + " name has not been initialized");
    }
    if (isotropic_exchange_parameters_names_.empty()) {
        isotropic_exchange_parameters_names_.resize(
            number_of_spins_,
            std::vector<SymbolName>(number_of_spins_));
    }
    if (!isotropic_exchange_parameters_names_[center_a][center_b].name.empty()) {
        throw std::invalid_argument("This parameter has been already specified");
    }

    if (symbols_[symbol_name].type_enum != SymbolTypeEnum::not_specified
        && symbols_[symbol_name].type_enum != SymbolTypeEnum::J) {
        throw std::invalid_argument(
            symbol_name.name + " has been already specified as not J parameter");
    }
    if (symbols_[symbol_name].type_enum == SymbolTypeEnum::not_specified) {
        symbols_[symbol_name].type_enum = SymbolTypeEnum::J;
    }

    isotropic_exchange_parameters_names_[center_a][center_b] = symbol_name;
    isotropic_exchange_parameters_names_[center_b][center_a] = symbol_name;
}

bool Symbols::hasIsotropicExchangeParameters() const {
    return isotropic_exchange_parameters_values_ != nullptr;
}

const Symbols::SymbolName Symbols::addSymbol(
    const std::string& name_string,
    double initial_value,
    bool is_changeable,
    SymbolTypeEnum type_enum) {
    SymbolName symbol_name = {name_string};
    if (symbols_.find(symbol_name) != symbols_.end()) {
        throw std::invalid_argument(symbol_name.name + " name has been already initialized");
    }
    symbols_[symbol_name] = {initial_value, is_changeable, type_enum};
    return symbol_name;
}

void Symbols::addGFactor(const SymbolName& symbol_name, size_t center_a) {
    if (symbols_.find(symbol_name) == symbols_.end()) {
        throw std::invalid_argument(symbol_name.name + " name has not been initialized");
    }

    if (symbols_[symbol_name].type_enum != SymbolTypeEnum::not_specified
        && symbols_[symbol_name].type_enum != SymbolTypeEnum::g_factor) {
        throw std::invalid_argument(
            symbol_name.name + " has been already specified as not g factor parameter");
    }
    if (symbols_[symbol_name].type_enum == SymbolTypeEnum::not_specified) {
        symbols_[symbol_name].type_enum = SymbolTypeEnum::g_factor;
    }

    if (g_factor_names_.empty()) {
        g_factor_names_.resize(number_of_spins_);
    }

    g_factor_names_[center_a] = symbol_name;
}

std::shared_ptr<const DenseVector> Symbols::constructGFactorParameters() {
    if (g_factor_names_.empty()) {
        throw std::length_error("g factor parameters have not been initialized");
    }

    for (const auto& symbol_name : g_factor_names_) {
        if (symbol_name.name.empty()) {
            throw std::invalid_argument("Not all g factors were initialized");
        }
    }

    if (g_factor_values_ == nullptr) {
        g_factor_values_ = std::make_shared<DenseVector>();
        g_factor_values_->resize(number_of_spins_);
    }

    updateGFactorParameters();

    return g_factor_values_;
}

bool Symbols::symmetry_consistence(const group::Group& group) const {
    // isotropic exchange parameters:
    if (!isotropic_exchange_parameters_names_.empty()) {
        for (const auto& element : group.elements_) {
            std::vector<std::vector<SymbolName>> permutated_symbols(
                element.size(),
                std::vector<SymbolName>(element.size()));
            for (size_t i = 0; i < element.size(); ++i) {
                for (size_t j = 0; j < element.size(); ++j) {
                    permutated_symbols[element[i]][element[j]] =
                        isotropic_exchange_parameters_names_[i][j];
                }
            }
            if (permutated_symbols != isotropic_exchange_parameters_names_) {
                return false;
            }
        }
    }
    if (!g_factor_names_.empty()) {
        for (const auto& element : group.elements_) {
            std::vector<SymbolName> permutated_symbols(element.size());
            for (size_t i = 0; i < element.size(); ++i) {
                permutated_symbols[element[i]] = g_factor_names_[i];
            }
            if (permutated_symbols != g_factor_names_) {
                return false;
            }
        }
    }
    return true;
}

std::shared_ptr<const DenseMatrix>
Symbols::constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name) {
    if (symbols_.find(symbol_name) == symbols_.end()) {
        throw std::invalid_argument(symbol_name.name + " name has not been initialized");
    }

    if (symbols_[symbol_name].type_enum != SymbolTypeEnum::J) {
        throw std::invalid_argument(symbol_name.name + " has been specified as not J parameter");
    }

    if (isotropic_exchange_parameters_names_.empty()) {
        throw std::length_error("Isotropic exchange parameters has not been initialized");
    }

    auto ptr_to_derivative = std::make_shared<DenseMatrix>();
    ptr_to_derivative->resize_with_nans(number_of_spins_, number_of_spins_);

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            if (isotropic_exchange_parameters_names_[i][j] != symbol_name) {
                value = NAN;
            } else {
                value = 1;
            }
            ptr_to_derivative->assign_to_position(value, i, j);
        }
    }

    isotropic_exchange_derivatives_values_.push_back(ptr_to_derivative);

    return ptr_to_derivative;
}

std::vector<Symbols::SymbolName> Symbols::getChangeableNames(SymbolTypeEnum type_enum) const {
    std::vector<SymbolName> answer;
    for (const auto& [symbol_name, symbol_data] : symbols_) {
        if (symbol_data.is_changeable && symbol_data.type_enum == type_enum) {
            answer.push_back(symbol_name);
        }
    }
    return answer;
}

bool Symbols::isAllGFactorsEqual() const {
    if (g_factor_names_.empty()) {
        throw std::length_error("g factor parameters have not been initialized");
    }
    for (size_t i = 1; i < g_factor_names_.size(); ++i) {
        if (g_factor_names_[i] != g_factor_names_[0]) {
            return false;
        }
    }
    return true;
}

std::vector<Symbols::SymbolName> Symbols::getChangeableNames() const {
    std::vector<SymbolName> answer;
    for (const auto& [symbol_name, symbol_data] : symbols_) {
        if (symbol_data.is_changeable) {
            answer.push_back(symbol_name);
        }
    }
    return answer;
}

double Symbols::getValueOfName(const SymbolName& symbol_name) const {
    return symbols_.at(symbol_name).value;
}

void Symbols::setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value) {
    SymbolData symbol_data = symbols_[symbol_name];
    if (!symbol_data.is_changeable) {
        throw std::invalid_argument("Cannot change value of unchangeable symbol");
    }
    symbol_data.value = new_value;
    symbols_[symbol_name] = symbol_data;
}

void Symbols::updateIsotropicExchangeParameters() {
    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            if (isotropic_exchange_parameters_names_[i][j].name.empty()) {
                value = NAN;
            } else {
                value = symbols_[isotropic_exchange_parameters_names_[i][j]].value;
            }
            isotropic_exchange_parameters_values_->assign_to_position(value, i, j);
        }
    }
}

void Symbols::updateGFactorParameters() {
    for (size_t i = 0; i < number_of_spins_; ++i) {
        double value = symbols_[g_factor_names_[i]].value;
        g_factor_values_->assign_to_position(value, i);
    }
}

bool Symbols::SymbolName::operator<(const Symbols::SymbolName& rhs) const {
    return name < rhs.name;
}
bool Symbols::SymbolName::operator>(const Symbols::SymbolName& rhs) const {
    return rhs < *this;
}
bool Symbols::SymbolName::operator<=(const Symbols::SymbolName& rhs) const {
    return !(rhs < *this);
}
bool Symbols::SymbolName::operator>=(const Symbols::SymbolName& rhs) const {
    return !(*this < rhs);
}
bool Symbols::SymbolName::operator==(const Symbols::SymbolName& rhs) const {
    return name == rhs.name;
}
bool Symbols::SymbolName::operator!=(const Symbols::SymbolName& rhs) const {
    return !(rhs == *this);
}
}  // namespace symbols