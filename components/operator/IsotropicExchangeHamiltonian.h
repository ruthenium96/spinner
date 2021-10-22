#ifndef JULY_ISOTROPICEXCHANGEHAMILTONIAN_H
#define JULY_ISOTROPICEXCHANGEHAMILTONIAN_H

#include "Interaction.h"

class IsotropicExchangeHamiltonian : public TwoCenterTerm {
public:

    explicit IsotropicExchangeHamiltonian(arma::dmat isotropic_exchange_parameters);

    void construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                   uint32_t index_of_vector, uint32_t center_a, uint32_t center_b) const override;

    [[nodiscard]] arma::dmat get_parameters() const override;

private:
    // TODO: abstract from arma::dmat
    arma::dmat isotropic_exchange_parameters;
};


#endif //JULY_ISOTROPICEXCHANGEHAMILTONIAN_H
