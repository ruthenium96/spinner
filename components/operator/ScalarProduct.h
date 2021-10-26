#ifndef JULY_SCALARPRODUCT_H
#define JULY_SCALARPRODUCT_H

#include "components/operator/Interaction.h"

class ScalarProduct : public TwoCenterTerm {
public:

    explicit ScalarProduct(size_t number_of_spins);
    explicit ScalarProduct(arma::dmat isotropic_exchange_parameters);

    void construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                   uint32_t index_of_vector, uint32_t center_a, uint32_t center_b) const override;

    [[nodiscard]] arma::dmat get_parameters() const override;

private:
    // TODO: abstract from arma::dmat
    arma::dmat coefficients;

    void scalar_product_total(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                              uint32_t index_of_vector, uint32_t center_a, uint32_t center_b, double factor) const;

    void scalar_product_nondiagonal_part(LexicograficalMatrix& matrix_in_lexicografical_basis,
                                         const spaces::LexicographicIndexConverter& converter,
                                         uint32_t index_of_vector, uint32_t plus_center, uint32_t minus_center,
                                         uint32_t projection_of_plus_center, uint32_t projection_of_minus_center,
                                         double factor) const;
};


#endif //JULY_SCALARPRODUCT_H
