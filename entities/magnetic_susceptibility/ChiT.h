#ifndef JULY_CHIT_H
#define JULY_CHIT_H

#include "EnsembleAverager.h"
#include "entities/data_structures/DenseMatrix.h"

// TODO: is it true?
constexpr double bohr_magneton = 0.67171388;

namespace magnetic_susceptibility {

enum QuantityEnum {
    mu_squared_in_bohr_magnetons_squared,
    chiT_in_cm_cubed_kelvin_per_mol,
};

struct MuSquaredData {
    double temperature;
    double mu_squared;
};

// mu^2 = mu_B^2 * g_{iso}^2 * <S^2_{total}>
//class SSquared {
//  public:
//    CalculatedCurve(double g_factor, DenseVector energy, DenseVector degeneracy, DenseVector s_squared);
//
//  private:
//    double g_factor_;
//    DenseVector energy_;
//    DenseVector s_squared_;
//    DenseVector degeneracy_;
//};

class ExperimentalCurve {
  private:
    QuantityEnum initial_quantity_type;
};

}  // namespace magnetic_susceptibility
#endif  //JULY_CHIT_H
