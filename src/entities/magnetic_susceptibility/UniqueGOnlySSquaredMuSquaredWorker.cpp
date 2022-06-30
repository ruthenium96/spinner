#include "UniqueGOnlySSquaredMuSquaredWorker.h"
namespace magnetic_susceptibility {

UniqueGOnlySSquaredMuSquaredWorker::UniqueGOnlySSquaredMuSquaredWorker(
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& energy,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& degeneracy,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& s_squared,
    double g_unique) :
    MuSquaredWorker(std::move(energy), std::move(degeneracy)),
    s_squared_(std::move(s_squared)),
    g_unique_(g_unique) {}

double UniqueGOnlySSquaredMuSquaredWorker::calculateTheoreticalMuSquared(double temperature) const {
    double s_squared_averaged = ensemble_averager_.ensemble_average(s_squared_, temperature);
    return g_unique_ * g_unique_ * s_squared_averaged;
}

std::vector<ValueAtTemperature> UniqueGOnlySSquaredMuSquaredWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    std::unique_ptr<quantum::linear_algebra::AbstractVector>&& derivative_value) const {
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> derivatives(temperatures.size());
    if (symbol_type == model::symbols::SymbolTypeEnum::g_factor) {
        throw std::invalid_argument("g-factor should be passed without Densevector.");
    } else {
        // d(mu_squared)/da = d(g^2*<S^2>)/da = g^2*d(<S^2>)/dg = g^2*(<S^2>*<dE/da>-<S^2*dE/da>)/T
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double first_term = ensemble_averager_.ensemble_average(s_squared_, temperatures[i])
                * ensemble_averager_.ensemble_average(derivative_value, temperatures[i]);
            double second_term = ensemble_averager_.ensemble_average(
                s_squared_->element_wise_multiplication(derivative_value),
                temperatures[i]);
            double value = g_unique_ * g_unique_ * (first_term - second_term) / temperatures[i];
            derivatives[i] = {temperatures[i], value};
        }
    }
    return derivatives;
}

std::vector<ValueAtTemperature> UniqueGOnlySSquaredMuSquaredWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type) const {
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> derivatives(temperatures.size());
    if (symbol_type == model::symbols::SymbolTypeEnum::g_factor) {
        // d(mu_squared)/dg = d(g^2*<S^2>)/dg = d(g^2)/dg*<S^2> + g^2*d(<S^2>)/dg = 2g*<S^2>
        for (size_t i = 0; i < temperatures.size(); ++i) {
            double value =
                2 * g_unique_ * ensemble_averager_.ensemble_average(s_squared_, temperatures[i]);
            derivatives[i] = {temperatures[i], value};
        }
    } else {
        throw std::invalid_argument("Only g-factor can be passed without Densevector.");
    }
    return derivatives;
}

}  // namespace magnetic_susceptibility
