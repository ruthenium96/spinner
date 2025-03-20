#ifndef SPINNER_ABSTRACTWORKER_H
#define SPINNER_ABSTRACTWORKER_H

#include "src/common/Quantity.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/model/symbols/SymbolName.h"
#include "src/model/symbols/SymbolicWorker.h"

namespace magnetic_susceptibility::worker {
class AbstractWorker {
  public:
    // Different cases lead to different approaches of calculating <S2>.
    virtual double calculateTheoreticalMuSquared(double temperature) const = 0;
    virtual std::shared_ptr<ExperimentalValuesWorker> getExperimentalValuesWorker() = 0;
    virtual std::shared_ptr<const ExperimentalValuesWorker> getExperimentalValuesWorker() const = 0;
    virtual void setExperimentalValuesWorker(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker) = 0;

    // d<A>/dsymbol
    virtual std::vector<ValueAtTemperature> calculateDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        model::symbols::SymbolName symbol_name) const = 0;            

    virtual ~AbstractWorker() = default;
};
}  // namespace magnetic_susceptibility::worker

#endif  //SPINNER_ABSTRACTWORKER_H
