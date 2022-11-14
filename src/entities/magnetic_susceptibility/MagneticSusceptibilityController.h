#ifndef SPINNER_MAGNETICSUSCEPTIBILITYCONTROLLER_H
#define SPINNER_MAGNETICSUSCEPTIBILITYCONTROLLER_H

#include "src/entities/data_structures/AbstractVector.h"
#include "src/entities/magnetic_susceptibility/assistant/ExperimentalValuesWorker.h"
#include "src/model/symbols/SymbolicWorker.h"
#include "worker/AbstractWorker.h"

// TODO: is it true?
//constexpr double bohr_magneton = 0.67171388;

namespace magnetic_susceptibility {

// Class for calculating mu^2 and dmu^2/dparameter.
// Different types of Hamiltonian and expressions in worker::AbstractWorker.
class MagneticSusceptibilityController {
  public:
    explicit MagneticSusceptibilityController(std::unique_ptr<worker::AbstractWorker>&& worker);
    void initializeExperimentalValues(
        const std::shared_ptr<ExperimentalValuesWorker>& experimental_values_worker);
    std::vector<ValueAtTemperature> getTheoreticalValues() const;
    double calculateTheoreticalMuSquared(double temperature) const;
    double calculateResidualError() const;
    // These functions just call suitable virtual function calculateDerivative.
    // Calculates dR^2/dsymbol.
    double calculateTotalDerivative(
        model::symbols::SymbolTypeEnum symbol_type,
        std::map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractVector>>
            values_derivatives_map) const;

  private:
    // dot product with some checks:
    double multiplyExperimentalAndTheoreticalDerivatives(
        std::vector<ValueAtTemperature> theoretical_derivative) const;
    std::unique_ptr<worker::AbstractWorker> worker_;
};

}  // namespace magnetic_susceptibility
#endif  //SPINNER_MAGNETICSUSCEPTIBILITYCONTROLLER_H
