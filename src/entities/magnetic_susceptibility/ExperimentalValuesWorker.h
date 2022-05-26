#ifndef SPINNER_EXPERIMENTALVALUESWORKER_H
#define SPINNER_EXPERIMENTALVALUESWORKER_H

#include <vector>

namespace magnetic_susceptibility {

constexpr double mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol =
    0.125048612;

struct ValueAtTemperature {
    double temperature;
    double value;
};

enum ExperimentalValuesEnum {
    mu_in_bohr_magnetons,
    mu_squared_in_bohr_magnetons_squared,
    chiT_in_cm_cubed_kelvin_per_mol,
};

class ExperimentalValuesWorker {
  public:
    // TODO: can we also add weight function?
    ExperimentalValuesWorker(
        std::vector<ValueAtTemperature> experimental_values,
        ExperimentalValuesEnum experimental_values_type,
        double number_of_centers_ratio);

    std::vector<double> getTemperatures() const;
    void setTheoreticalMuSquared(std::vector<ValueAtTemperature> theoretical_mu_squared);

    std::vector<ValueAtTemperature> getTheoreticalValues() const;
    std::vector<ValueAtTemperature> getExperimentalMuSquared() const;

    double calculateResidualError() const;
    std::vector<ValueAtTemperature> calculateDerivative() const;

  private:
    void throw_exception_if_theoretical_mu_squared_was_not_initialized() const;

    std::vector<ValueAtTemperature> experimental_mu_squared_;
    std::vector<ValueAtTemperature> theoretical_mu_squared_;
    ExperimentalValuesEnum experimental_values_type_;
    double number_of_centers_ratio_;
};
}  // namespace magnetic_susceptibility
#endif  //SPINNER_EXPERIMENTALVALUESWORKER_H
