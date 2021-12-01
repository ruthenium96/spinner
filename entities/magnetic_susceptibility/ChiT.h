#ifndef JULY_CHIT_H
#define JULY_CHIT_H

#include <optional>

#include "EnsembleAverager.h"
#include "ExperimentalValuesWorker.h"
#include "entities/data_structures/DenseMatrix.h"

// TODO: is it true?
//constexpr double bohr_magneton = 0.67171388;

namespace magnetic_susceptibility {

class MuSquaredWorker {
  public:
    MuSquaredWorker(DenseVector&& energy, DenseVector&& degeneracy);

    void initializeExperimentalValues(
        std::vector<ValueAtTemperature> experimental_values,
        ExperimentalValuesEnum experimental_quantity_type,
        double number_of_centers_ratio);

    void calculateTheoreticalValues();

  private:
    //    class MuSquaredCalculator {};
    std::optional<ExperimentalValuesWorker> experimental_values_worker_;
    const EnsembleAverager ensemble_averager_;
    //    std::vector<ValueAtTemperature> theoretical_mu_squared;
};

// mu^2 = mu_B^2 * g_{iso}^2 * <S^2_{total}>

}  // namespace magnetic_susceptibility
#endif  //JULY_CHIT_H
