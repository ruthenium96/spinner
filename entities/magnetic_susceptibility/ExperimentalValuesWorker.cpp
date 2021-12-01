#include "ExperimentalValuesWorker.h"

namespace magnetic_susceptibility {

ExperimentalValuesWorker::ExperimentalValuesWorker(
    std::vector<ValueAtTemperature> experimental_values,
    ExperimentalValuesEnum experimental_values_type,
    double number_of_centers_ratio) {
    experimental_values_type_ = experimental_values_type;

    if (experimental_values_type == mu_squared_in_bohr_magnetons_squared) {
        // DO NOTHING
    } else if (experimental_values_type == mu_in_bohr_magnetons) {
        for (ValueAtTemperature& v : experimental_values) {
            v.value = v.value * v.value;
        }
    } else if (experimental_values_type == chiT_in_cm_cubed_kelvin_per_mol) {
        for (ValueAtTemperature& v : experimental_values) {
            v.value /= mu_squared_in_bohr_magnetons_squared_to_chiT_in_cm_cubed_kelvin_per_mol;
        }
    }

    for (ValueAtTemperature& v : experimental_values) {
        v.value *= number_of_centers_ratio;
    }

    experimental_mu_squared = std::move(experimental_values);
}
}  // namespace magnetic_susceptibility