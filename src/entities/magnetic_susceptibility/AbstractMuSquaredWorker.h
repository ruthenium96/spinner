#ifndef SPINNER_ABSTRACTMUSQUAREDWORKER_H
#define SPINNER_ABSTRACTMUSQUAREDWORKER_H

#include "src/entities/data_structures/AbstractVector.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/model/symbols/Symbols.h"

// TODO: is it true?
//constexpr double bohr_magneton = 0.67171388;

namespace magnetic_susceptibility {

// Abstract class for calculating mu^2 and dmu^2/dparameter.
// Concrete classes correspond to different types of Hamiltonian or expressions.
class AbstractMuSquaredWorker {
  public:
    // Different cases lead to different approaches of calculating <S2>.
    // TODO: is <S2> required in all cases?
    virtual double calculateTheoreticalMuSquared(double temperature) const = 0;
    virtual std::shared_ptr<ExperimentalValuesWorker> getExperimentalValuesWorker() = 0;
    virtual std::shared_ptr<const ExperimentalValuesWorker> getExperimentalValuesWorker() const = 0;
    virtual void setExperimentalValuesWorker(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) = 0;

    // Some values do not change <A> values, so d<A>/dvalue = 0. Function for this case:
    virtual std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const = 0;
    // Some values change <A> values, so we need d<A>/dvalue. Function for this case:
    virtual std::vector<ValueAtTemperature>
    calculateDerivative(model::symbols::SymbolTypeEnum symbol_type) const = 0;

  public:
    void initializeExperimentalValues(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker);
    std::vector<ValueAtTemperature> getTheoreticalValues() const;
    double calculateResidualError() const;
    // These functions just call suitable virtual function calculateDerivative.
    // Some values change <A> values, so we need d<A>/dvalue. Function for this case:
    double calculateTotalDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const;
    // Some values do not change <A> values, so d<A>/dvalue = 0. Function for this case:
    double calculateTotalDerivative(model::symbols::SymbolTypeEnum symbol_type) const;
  protected:
    // dot product with some checks:
    double multiplyExperimentalAndTheoreticalDerivatives(
        std::vector<ValueAtTemperature> theoretical_derivative) const;
};

}  // namespace magnetic_susceptibility
#endif  //SPINNER_ABSTRACTMUSQUAREDWORKER_H
