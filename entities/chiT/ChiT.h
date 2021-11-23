#ifndef JULY_CHIT_H
#define JULY_CHIT_H

#include "entities/data_structures/DenseMatrix.h"

// TODO: is it true?
constexpr double bohr_magneton = 0.67171388;

namespace magnetic_susceptibility {

class ChiT {
  public:
    double at_temperature(double temperature);
    // chiT = mu_B^2 * g_{iso}^2 / 3 * <S^2_{total}>
    ChiT(double g_factor, DenseVector energy, DenseVector degeneracy, DenseVector s_squared);

  private:
    double ensemble_average(
        double temperature,
        const DenseVector& energy,
        const DenseVector& degeneracy,
        const DenseVector& value);
    double g_factor_;
    DenseVector energy_;
    DenseVector s_squared_;
    DenseVector degeneracy_;
};

}  // namespace magnetic_susceptibility
#endif  //JULY_CHIT_H
