#include "UniqueGWorker.h"
namespace magnetic_susceptibility::worker {

UniqueGWorker::UniqueGWorker(
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& energy,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& degeneracy,
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>&& quantity,
    double g_unique) :
    BasicWorker(std::move(energy), std::move(degeneracy)),
    quantity_(std::move(quantity)),
    g_unique_(g_unique) {}

double UniqueGWorker::calculateTheoreticalMuSquared(double temperature) const {
    double s_squared_averaged = ensemble_averager_.ensemble_average(quantity_, temperature);
    return g_unique_ * g_unique_ * s_squared_averaged;
}

std::vector<ValueAtTemperature> UniqueGWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::map<common::QuantityEnum, std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>
        values_derivatives_map) const {
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> derivatives(temperatures.size());
    if (symbol_type == model::symbols::SymbolTypeEnum::g_factor) {
        // d(mu_squared)/dg = d(g^2*<S^2>)/dg = d(g^2)/dg*<S^2> + g^2*d(<S^2>)/dg = 2g*<S^2>
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double value =
                2 * g_unique_ * ensemble_averager_.ensemble_average(quantity_, temperatures[i]);
            derivatives[i] = {temperatures[i], value};
        }
    } else {
        // d(mu_squared)/da = d(g^2*<S^2>)/da = g^2*d(<S^2>)/da = g^2*(<S^2>*<dE/da>-<S^2*dE/da>)/T
        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& energy_derivative =
            values_derivatives_map[common::Energy];
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double first_term = ensemble_averager_.ensemble_average(quantity_, temperatures[i])
                * ensemble_averager_.ensemble_average(energy_derivative, temperatures[i]);
            double second_term = ensemble_averager_.ensemble_average(
                quantity_->element_wise_multiplication(energy_derivative),
                temperatures[i]);
            double value = g_unique_ * g_unique_ * (first_term - second_term) / temperatures[i];
            derivatives[i] = {temperatures[i], value};
        }
    }
    return derivatives;
}

}  // namespace magnetic_susceptibility
