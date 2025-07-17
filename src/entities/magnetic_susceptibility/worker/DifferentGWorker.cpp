#include "DifferentGWorker.h"
#include "src/common/Quantity.h"
#include "src/common/UncertainValue.h"
#include "src/model/symbols/SymbolName.h"

namespace magnetic_susceptibility::worker {
DifferentGWorker::DifferentGWorker(
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra,
    common::QuantityEnum quantity_enum_for_averaging, 
    double quantity_factor) :
    BasicWorker(flattenedSpectra),
    flattenedSpectra_(flattenedSpectra),
    quantity_enum_for_averaging_(quantity_enum_for_averaging),
    quantity_factor_(quantity_factor) {}

common::UncertainValue DifferentGWorker::calculateTheoreticalMuSquared(double temperature) const {
    auto quantity = flattenedSpectra_->getFlattenSpectrum(quantity_enum_for_averaging_).value();
    auto quantity_averaged =
        quantity_factor_ * ensemble_averager_->ensemble_average(quantity, temperature);
    return quantity_averaged;
}

std::vector<ValueAtTemperature> DifferentGWorker::calculateDerivative(
    model::symbols::SymbolTypeEnum symbol_type,
    model::symbols::SymbolName symbol_name) const {
    auto quantity = flattenedSpectra_->getFlattenSpectrum(quantity_enum_for_averaging_).value();
    std::vector<double> temperatures = experimental_values_worker_.value()->getTemperatures();
    std::vector<ValueAtTemperature> derivatives(temperatures.size());
    // let A = (\sum_i g_i S_{iz})^2
    if (symbol_type == model::symbols::SymbolTypeEnum::g_factor) {
        // if energy does not depend on g factors:
        // d(mu_squared)/dg = d<A>/dg = <dA/dg>
        auto quantity_derivative = 
            flattenedSpectra_->getFlattenDerivativeSpectrum(quantity_enum_for_averaging_, symbol_name).value();
        for (size_t i = 0; i < temperatures.size(); ++i) {
            auto value =
                quantity_factor_ * ensemble_averager_->ensemble_average(quantity_derivative, temperatures[i]);
            derivatives[i] = {temperatures[i], value};
        }
    } else {
        // d(mu_squared)/da = d(<A>)/da = (<A>*<dE/da>-<A*dE/da>)/T
        auto energy_derivative =
            flattenedSpectra_->getFlattenDerivativeSpectrum(common::Energy, symbol_name).value();
        for (size_t i = 0; i < temperatures.size(); ++i) {
            auto first_term = ensemble_averager_->ensemble_average(quantity, temperatures[i])
                * ensemble_averager_->ensemble_average(energy_derivative, temperatures[i]);
            auto quantity_derivative_product = flattenedSpectra_->getFlattenDerivativeProductSpectrum(quantity_enum_for_averaging_, common::Energy, symbol_name).value();
            auto second_term = ensemble_averager_->ensemble_average(quantity_derivative_product, temperatures[i]);
            auto value = quantity_factor_ * (first_term - second_term) / temperatures[i];
            derivatives[i] = {temperatures[i], value};
        }
    }
    return derivatives;
}

}  // namespace magnetic_susceptibility::worker