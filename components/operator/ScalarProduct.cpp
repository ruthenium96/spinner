#include "ScalarProduct.h"


ScalarProduct::ScalarProduct(size_t number_of_spins) {
    coefficients.resize(number_of_spins, number_of_spins);
    for (size_t i = 0; i < number_of_spins; ++i) {
        for (size_t j = i + 1; j < number_of_spins; ++j) {
            // TODO: it does not looks good
            coefficients(i, j) = -1;
            coefficients(j, i) = -1;
        }
    }
}

ScalarProduct::ScalarProduct(arma::dmat parameters) : coefficients(std::move(parameters)) {
}

void ScalarProduct::construct(LexicographicSparseMatrix& matrix_in_lexicografical_basis, uint32_t index_of_vector, uint32_t center_a, uint32_t center_b) const {
    if (!std::isnan(coefficients(center_a, center_b))) {
        double factor = -2 * coefficients(center_a, center_b);
        matrix_in_lexicografical_basis.add_scalar_product(index_of_vector, center_a, center_b, factor);
    }
}

arma::dmat ScalarProduct::get_parameters() const {
    return coefficients;
}
