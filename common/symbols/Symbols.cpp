#include "Symbols.h"

#include <cmath>
#include <utility>

namespace symbols {
Symbols::Symbols(size_t number_of_spins) : number_of_spins_(number_of_spins) {}

std::shared_ptr<const DenseMatrix> symbols::Symbols::constructIsotropicExchangeParameters() {
    if (isotropic_exchange_parameters_values_ == nullptr) {
        isotropic_exchange_parameters_values_ = std::make_shared<DenseMatrix>();
        isotropic_exchange_parameters_values_->resize_with_nans(number_of_spins_, number_of_spins_);
    }

    for (size_t i = 0; i < number_of_spins_; ++i) {
        for (size_t j = 0; j < number_of_spins_; ++j) {
            double value;
            if (isotropic_exchange_parameters_symbols_[i][j].empty()) {
                value = NAN;
            } else {
                value = name_to_value_[isotropic_exchange_parameters_symbols_[i][j]];
            }
            isotropic_exchange_parameters_values_->assign_to_position(value, i, j);
        }
    }

    return isotropic_exchange_parameters_values_;
}

void Symbols::addIsotropicExchange(
    const std::string& symbol_name,
    size_t center_a,
    size_t center_b) {
    if (isotropic_exchange_parameters_symbols_.empty()) {
        isotropic_exchange_parameters_symbols_.resize(
            number_of_spins_,
            std::vector<std::string>(number_of_spins_));
    }
    if (!isotropic_exchange_parameters_symbols_[center_a][center_b].empty()) {
        throw std::invalid_argument("This parameter has been already specified");
    }

    isotropic_exchange_parameters_symbols_[center_a][center_b] = symbol_name;
    isotropic_exchange_parameters_symbols_[center_b][center_a] = symbol_name;
}

bool Symbols::hasIsotropicExchangeParameters() const {
    return isotropic_exchange_parameters_values_ != nullptr;
}

void Symbols::addSymbol(const std::string& name, double initial_value, bool is_changeable) {
    name_to_value_[name] = initial_value;
    if (is_changeable) {
        changeable_symbols_.emplace_back(name);
    }
}

}  // namespace symbols