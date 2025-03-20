#include "DifferentGWorker.h"
#include "src/common/Quantity.h"

namespace magnetic_susceptibility::worker {
DifferentGWorker::DifferentGWorker(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra,
    common::QuantityEnum quantity_enum_for_averaging, 
    double quantity_factor) :
    BasicWorker(flattenedSpectra),
    flattenedSpectra_(flattenedSpectra),
    quantity_enum_for_averaging_(quantity_enum_for_averaging),
    quantity_factor_(quantity_factor) {}

double DifferentGWorker::calculateTheoreticalMuSquared(double temperature) const {
    const auto& quantity = flattenedSpectra_->getFlattenSpectrum(quantity_enum_for_averaging_).value().get();
    double quantity_averaged =
        quantity_factor_ * ensemble_averager_.ensemble_average(quantity, temperature);
    return quantity_averaged;
}

std::vector<ValueAtTemperature> DifferentGWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
        values_derivatives_map) const {
    const auto& quantity = flattenedSpectra_->getFlattenSpectrum(quantity_enum_for_averaging_).value().get();
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> derivatives(temperatures.size());
    // let A = (\sum_i g_i S_{iz})^2
    if (symbol_type == model::symbols::SymbolTypeEnum::g_factor) {
        // if energy does not depend on g factors:
        // d(mu_squared)/dg = d<A>/dg = <dA/dg>
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&
            g_sz_squared_derivative = values_derivatives_map[common::gSz_total_squared];
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double value =
                quantity_factor_ * ensemble_averager_.ensemble_average(g_sz_squared_derivative, temperatures[i]);
            derivatives[i] = {temperatures[i], value};
        }
    } else {
        // d(mu_squared)/da = d(<A>)/da = (<A>*<dE/da>-<A*dE/da>)/T
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& energy_derivative =
            values_derivatives_map[common::Energy];
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double first_term = ensemble_averager_.ensemble_average(quantity, temperatures[i])
                * ensemble_averager_.ensemble_average(energy_derivative, temperatures[i]);
            double second_term = ensemble_averager_.ensemble_average(
                quantity->element_wise_multiplication(energy_derivative),
                temperatures[i]);
            double value = quantity_factor_ * (first_term - second_term) / temperatures[i];
            derivatives[i] = {temperatures[i], value};
        }
    }
    return derivatives;
}

}  // namespace magnetic_susceptibility::worker