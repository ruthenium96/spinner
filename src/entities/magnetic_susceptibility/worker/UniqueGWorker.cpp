#include "UniqueGWorker.h"
#include "src/common/UncertainValue.h"
#include "src/model/symbols/SymbolName.h"

namespace magnetic_susceptibility::worker {

UniqueGWorker::UniqueGWorker(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra,
    std::function<double()> g_unique_getter,
    common::QuantityEnum quantity_enum_for_averaging, 
    double quantity_factor) :
    BasicWorker(flattenedSpectra),
    flattenedSpectra_(flattenedSpectra),
    g_unique_getter_(g_unique_getter),
    quantity_enum_for_averaging_(quantity_enum_for_averaging),
    quantity_factor_(quantity_factor) {}

common::UncertainValue UniqueGWorker::calculateTheoreticalMuSquared(double temperature) const {
    auto quantity = flattenedSpectra_->getFlattenSpectrum(quantity_enum_for_averaging_).value();
    auto quantity_averaged = ensemble_averager_->ensemble_average(quantity, temperature);
    return g_unique_getter_() * g_unique_getter_() * quantity_averaged * quantity_factor_;
}

std::vector<ValueAtTemperature> UniqueGWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    model::symbols::SymbolName symbol_name) const {
    auto quantity = flattenedSpectra_->getFlattenSpectrum(quantity_enum_for_averaging_).value();
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> derivatives(temperatures.size());
    if (symbol_type == model::symbols::SymbolTypeEnum::g_factor) {
        // d(mu_squared)/dg = d(g^2*<S^2>)/dg = d(g^2)/dg*<S^2> + g^2*d(<S^2>)/dg = 2g*<S^2>
        for (size_t i = 0; i < temperatures.size(); ++i) {
            auto quantity_averaged = ensemble_averager_->ensemble_average(quantity, temperatures[i]);
            auto value = 2 * g_unique_getter_() * quantity_averaged * quantity_factor_;
            derivatives[i] = {temperatures[i], value};
        }
    } else {
        // d(mu_squared)/da = d(g^2*<S^2>)/da = g^2*d(<S^2>)/da = g^2*(<S^2>*<dE/da>-<S^2*dE/da>)/T
        auto energy_derivative =
            flattenedSpectra_->getFlattenDerivativeSpectrum(common::Energy, symbol_name).value();
        for (size_t i = 0; i < temperatures.size(); ++i) {
            auto first_term = ensemble_averager_->ensemble_average(quantity, temperatures[i])
                * ensemble_averager_->ensemble_average(energy_derivative, temperatures[i]);
            auto quantity_derivative_product = flattenedSpectra_->getFlattenDerivativeProductSpectrum(quantity_enum_for_averaging_, common::Energy, symbol_name).value();
            auto second_term = ensemble_averager_->ensemble_average(quantity_derivative_product, temperatures[i]);
            auto value = g_unique_getter_() * g_unique_getter_() * quantity_factor_ * 
                (first_term - second_term) / temperatures[i];
            derivatives[i] = {temperatures[i], value};
        }
    }
    return derivatives;
}

}  // namespace magnetic_susceptibility
