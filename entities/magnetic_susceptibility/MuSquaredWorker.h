#ifndef JULY_MUSQUAREDWORKER_H
#define JULY_MUSQUAREDWORKER_H

#include <optional>

#include "EnsembleAverager.h"
#include "ExperimentalValuesWorker.h"
#include "common/symbols/Symbols.h"
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

    std::vector<ValueAtTemperature> getTheoreticalValues() const;

    virtual double theory_at_temperature(double temperature) const = 0;

    double calculateResidualError() const;
    double calculateTotalDerivative(
        symbols::SymbolTypeEnum symbol_type,
        DenseVector&& derivative_value) const;
    double calculateTotalDerivative(symbols::SymbolTypeEnum symbol_type) const;

  protected:
    const EnsembleAverager ensemble_averager_;
    virtual std::vector<ValueAtTemperature> calculateDerivative(
        symbols::SymbolTypeEnum symbol_type,
        DenseVector&& derivative_value) const = 0;
    virtual std::vector<ValueAtTemperature>
    calculateDerivative(symbols::SymbolTypeEnum symbol_type) const = 0;

    double multiplyExperimentalAndTheoreticalDerivatives(
        std::vector<ValueAtTemperature> theoretical_derivative) const;

    std::optional<ExperimentalValuesWorker> experimental_values_worker_;
};

}  // namespace magnetic_susceptibility
#endif  //JULY_MUSQUAREDWORKER_H
