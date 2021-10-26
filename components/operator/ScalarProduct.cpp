#include "ScalarProduct.h"


ScalarProduct::ScalarProduct(size_t number_of_spins) {
    coefficients.resize(number_of_spins, number_of_spins);
    for (int i = 0; i < number_of_spins; ++i) {
        for (int j = i + 1; j < number_of_spins; ++j) {
            // TODO: it does not looks good
            coefficients(i, j) = -1;
            coefficients(j, i) = -1;
        }
    }
}

ScalarProduct::ScalarProduct(arma::dmat parameters) : coefficients(std::move(parameters)) {
}

void ScalarProduct::construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                                             uint32_t index_of_vector, uint32_t center_a, uint32_t center_b) const {
    if (!std::isnan(coefficients(center_a, center_b))) {
        double factor = -2 * coefficients(center_a, center_b);
        scalar_product_total(matrix_in_lexicografical_basis, converter, index_of_vector, center_a, center_b, factor);
    }
}

arma::dmat ScalarProduct::get_parameters() const {
    return coefficients;
}

void ScalarProduct::scalar_product_total(LexicograficalMatrix& matrix_in_lexicografical_basis,
                          const spaces::LexicographicIndexConverter& converter,
                          uint32_t index_of_vector, uint32_t center_a,
                          uint32_t center_b, double factor) const {
    uint32_t projection_of_center_a = converter.convert_lex_index_to_one_sz_projection(index_of_vector, center_a);
    uint32_t projection_of_center_b = converter.convert_lex_index_to_one_sz_projection(index_of_vector, center_b);

    // (Sa, Sb) = Sax*Sbx + Say*Sby + Saz*Sbz = 0.5 * (Sa+Sb- + Sa-Sb+) + Saz*Sbz

    // Saz Sbz
    matrix_in_lexicografical_basis(index_of_vector, index_of_vector) += (projection_of_center_a - converter.get_spins()[center_a]) *
            (projection_of_center_b - converter.get_spins()[center_b]) * factor;
    // Sa+ Sb-
    scalar_product_nondiagonal_part(matrix_in_lexicografical_basis, converter, index_of_vector, center_a, center_b,
                                    projection_of_center_a, projection_of_center_b, factor);
    // Sa- Sb+
    scalar_product_nondiagonal_part(matrix_in_lexicografical_basis, converter, index_of_vector, center_b, center_a,
                                    projection_of_center_b, projection_of_center_a, factor);

}

void ScalarProduct::scalar_product_nondiagonal_part(LexicograficalMatrix& matrix_in_lexicografical_basis,
                                     const spaces::LexicographicIndexConverter& converter,
                                     uint32_t index_of_vector,
                                     uint32_t plus_center, uint32_t minus_center,
                                     uint32_t projection_of_plus_center, uint32_t projection_of_minus_center,
                                     double factor) const {
    if (projection_of_plus_center == converter.get_mults()[plus_center] - 1 || projection_of_minus_center == 0) {
        return;
    }
    unsigned long index_of_new_vector = converter.ladder_projection(
            converter.ladder_projection(
                    index_of_vector, plus_center, +1), minus_center, -1
                    );

    // projection m = number n - spin S
    // so S(S+1)-m(m+1) = (2S-n)(n+1)
    // so S(S+1)-m(m-1) = n(2S+1-n)
    double factor_a = (2 * converter.get_spins()[plus_center] - projection_of_plus_center) * (projection_of_plus_center + 1);
    double factor_b = projection_of_minus_center * (2 * converter.get_spins()[minus_center] + 1 - projection_of_minus_center);

    // TODO: fix plus-minus
    matrix_in_lexicografical_basis(index_of_vector, index_of_new_vector) += 0.5 * sqrt(factor_a * factor_b) * factor;
}