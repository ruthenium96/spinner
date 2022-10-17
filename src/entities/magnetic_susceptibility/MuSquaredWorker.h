#ifndef SPINNER_MUSQUAREDWORKER_H
#define SPINNER_MUSQUAREDWORKER_H

#include <optional>

#include "EnsembleAverager.h"
#include "ExperimentalValuesWorker.h"
#include "src/entities/data_structures/AbstractVector.h"
#include "src/model/symbols/Symbols.h"

// TODO: is it true?
//constexpr double bohr_magneton = 0.67171388;

namespace magnetic_susceptibility {

// Abstract class for calculating mu^2 and dmu^2/dparameter.
// Concrete classes correspond to different types of Hamiltonian.
class MuSquaredWorker {
  public:
    // Different cases lead to different approaches of calculating <S2>.
    // TODO: is <S2> required in all cases?
    virtual double calculateTheoreticalMuSquared(double temperature) const = 0;

  protected:
    // Some values do not change <A> values, so d<A>/dvalue = 0. Function for this case:
    virtual std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const = 0;
    // Some values change <A> values, so we need d<A>/dvalue. Function for this case:
    virtual std::vector<ValueAtTemperature>
    calculateDerivative(model::symbols::SymbolTypeEnum symbol_type) const = 0;

  public:
    MuSquaredWorker(
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy);

    void initializeExperimentalValues(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker);

    std::vector<ValueAtTemperature> getTheoreticalValues() const;

    double calculateResidualError() const;
    // These functions just call suitable virtual function calculateDerivative.
    double calculateTotalDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const;
    double calculateTotalDerivative(model::symbols::SymbolTypeEnum symbol_type) const;

  protected:
    const EnsembleAverager ensemble_averager_;

    // dot product with some checks:
    double multiplyExperimentalAndTheoreticalDerivatives(
        std::vector<ValueAtTemperature> theoretical_derivative) const;

    std::optional<std::shared_ptr<ExperimentalValuesWorker>> experimental_values_worker_;
};

}  // namespace magnetic_susceptibility
#endif  //SPINNER_MUSQUAREDWORKER_H
