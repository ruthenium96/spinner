#include "IsotropicExchangeHamiltonian.h"
#include "ScalarProduct.h"


IsotropicExchangeHamiltonian::IsotropicExchangeHamiltonian(arma::dmat parameters) : isotropic_exchange_parameters(std::move(parameters)) {
}

void IsotropicExchangeHamiltonian::construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                                             uint32_t index_of_vector, uint32_t center_a, uint32_t center_b) const {
    if (!std::isnan(isotropic_exchange_parameters(center_a, center_b))) {
        double factor = -2 * isotropic_exchange_parameters(center_a, center_b);
        scalar_product_total(matrix_in_lexicografical_basis, converter, index_of_vector, center_a, center_b, factor);
    }
}

arma::dmat IsotropicExchangeHamiltonian::get_parameters() const {
    return isotropic_exchange_parameters;
}