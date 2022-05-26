#ifndef JULY_MUSQUAREDWORKER_H
#define JULY_MUSQUAREDWORKER_H

#include <optional>

#include "EnsembleAverager.h"
#include "ExperimentalValuesWorker.h"
#include "src/entities/data_structures/DenseMatrix.h"
#include "src/model/symbols/Symbols.h"

// TODO: is it true?
//constexpr double bohr_magneton = 0.67171388;

namespace magnetic_susceptibility {

class MuSquaredWorker {
  public:
    MuSquaredWorker(DenseVector&& energy, DenseVector&& degeneracy);

    void initializeExperimentalValues(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker);

    std::vector<ValueAtTemperature> getTheoreticalValues() const;

    virtual double calculateTheoreticalMuSquared(double temperature) const = 0;

    double calculateResidualError() const;
    double calculateTotalDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        DenseVector&& derivative_value) const;
    double calculateTotalDerivative(model::symbols::SymbolTypeEnum symbol_type) const;

  protected:
    const EnsembleAverager ensemble_averager_;
    virtual std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        DenseVector&& derivative_value) const = 0;
    virtual std::vector<ValueAtTemperature>
    calculateDerivative(model::symbols::SymbolTypeEnum symbol_type) const = 0;

    double multiplyExperimentalAndTheoreticalDerivatives(
        std::vector<ValueAtTemperature> theoretical_derivative) const;

    std::optional<std::shared_ptr<ExperimentalValuesWorker>> experimental_values_worker_;
};

}  // namespace magnetic_susceptibility
#endif  //JULY_MUSQUAREDWORKER_H
