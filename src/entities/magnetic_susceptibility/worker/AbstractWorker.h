#ifndef SPINNER_ABSTRACTWORKER_H
#define SPINNER_ABSTRACTWORKER_H

#include "src/entities/data_structures/AbstractVector.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/model/symbols/Symbols.h"

namespace magnetic_susceptibility::worker {
class AbstractWorker {
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
};
}  // namespace magnetic_susceptibility::worker

#endif  //SPINNER_ABSTRACTWORKER_H
