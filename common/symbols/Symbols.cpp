#include "Symbols.h"

#include <cmath>
#include <utility>

namespace symbols {
Symbols::Symbols(size_t number_of_spins) : number_of_spins_(number_of_spins) {}

std::shared_ptr<const DenseMatrix> symbols::Symbols::getIsotropicExchangeParameters() const {
    return isotropic_exchange_parameters_;
}

void Symbols::addIsotropicExchange(double value, size_t center_a, size_t center_b) {
    if (center_b == center_a) {
        throw std::invalid_argument("Isotropic exchange takes place between different centers");
    }
    if (isotropic_exchange_parameters_ == nullptr) {
        isotropic_exchange_parameters_ = std::make_shared<DenseMatrix>();
        isotropic_exchange_parameters_->resize_with_nans(number_of_spins_, number_of_spins_);
    }

    if (!std::isnan(isotropic_exchange_parameters_->operator()(center_a, center_b))) {
        throw std::invalid_argument("This parameter has been already specified");
    }

    isotropic_exchange_parameters_->assign_to_position(value, center_a, center_b);
    isotropic_exchange_parameters_->assign_to_position(value, center_b, center_a);
}

bool Symbols::hasIsotropicExchangeParameters() const {
    return isotropic_exchange_parameters_ != nullptr;
}

}  // namespace symbols