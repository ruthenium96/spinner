#ifndef JULY_SCALARPRODUCT_H
#define JULY_SCALARPRODUCT_H

    // TODO: This function is also required for S^2 matrix building
void scalar_product_total(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                          uint32_t index_of_vector, uint32_t center_a, uint32_t center_b, double factor);

void scalar_product_nondiagonal_part(LexicograficalMatrix& matrix_in_lexicografical_basis,
                                     const spaces::LexicographicIndexConverter& converter,
                                     uint32_t index_of_vector, uint32_t plus_center, uint32_t minus_center,
                                     uint32_t projection_of_plus_center, uint32_t projection_of_minus_center,
                                     double factor);

void scalar_product_total(LexicograficalMatrix& matrix_in_lexicografical_basis,
                                                        const spaces::LexicographicIndexConverter& converter,
                                                        uint32_t index_of_vector, uint32_t center_a,
                                                        uint32_t center_b, double factor) {
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

void scalar_product_nondiagonal_part(LexicograficalMatrix& matrix_in_lexicografical_basis,
                                                                   const spaces::LexicographicIndexConverter& converter,
                                                                   uint32_t index_of_vector,
                                                                   uint32_t plus_center, uint32_t minus_center,
                                                                   uint32_t projection_of_plus_center, uint32_t projection_of_minus_center,
                                                                   double factor) {
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

#endif //JULY_SCALARPRODUCT_H
