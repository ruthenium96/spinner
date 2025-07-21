#include "FlattenedSpectra.h"

#include <functional>
#include <memory>
#include <optional>

#include "magic_enum.hpp"
#include "src/common/Quantity.h"
#include "src/eigendecompositor/AllQuantitiesGetter.h"
#include "src/entities/spectrum/Spectrum.h"

namespace {

OneOrMany<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> flatten(
    OneOrMany<SpectrumRef> spectrum,
    const quantum::linear_algebra::FactoriesList& factories) {
    return std::move(transform_one_or_many(
        std::function([factories](SpectrumRef spectrum){
            auto vector = factories.createVector();
            for (const auto& subspectrum_ref : spectrum.blocks) {
                const auto& subspectrum = subspectrum_ref.get();
                vector->concatenate_with(subspectrum.raw_data);
            }
            return std::move(vector);
        }), 
        spectrum));
    }
} // namespace

namespace eigendecompositor {

void FlattenedSpectra::updateValues(const AllQuantitiesGetter& allQuantitiesGetter,
    const quantum::linear_algebra::FactoriesList& factories) {
    flattenedSpectra_.clear();
    degeneracyValues_ = factories.createVector();
    for (const auto& quantity_enum : magic_enum::enum_values<common::QuantityEnum>()) {
        auto maybe_spectrum = allQuantitiesGetter.getSpectrum(quantity_enum);
        if (maybe_spectrum.has_value()) {
            auto spectrum = maybe_spectrum.value();
            flattenedSpectra_[quantity_enum] = flatten(spectrum, factories);
        }
    }
    flattenedWeights_ = transform_one_or_many(
        std::function([factories](std::vector<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> weights_all){
            auto vector = factories.createVector();
            for (auto ref_weights_of_block : weights_all) {
                const auto& weights_of_block = ref_weights_of_block.get();
                vector->concatenate_with(weights_of_block);
            }
            return std::move(vector);
        }),
        allQuantitiesGetter.getWeightsOfAllStates());
    apply_to_one_or_many(
        std::function([](const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& energy_spectrum){
            energy_spectrum.get()->subtract_minimum();
        }), 
        flattenedSpectra_[common::Energy]);
}

void FlattenedSpectra::updateDerivativeValues(const AllQuantitiesGetter& allQuantitiesGetter,
    const std::vector<model::symbols::SymbolName>& symbol_names,
    const quantum::linear_algebra::FactoriesList& factories) {
    flattenedDerivativeSpectra_.clear();
    for (const auto& quantity_enum : magic_enum::enum_values<common::QuantityEnum>()) {
        for (const auto& symbol_name : symbol_names) {
            auto maybe_spectrum = allQuantitiesGetter.getSpectrumDerivative(quantity_enum, symbol_name);
            if (maybe_spectrum.has_value()) {
                const auto& spectrum = maybe_spectrum.value();
                flattenedDerivativeSpectra_[{quantity_enum, symbol_name}] = flatten(spectrum, factories);
            }
        }
    }   
    for (const auto& [quantity_enum, spectrum] : flattenedSpectra_) {
        if (quantity_enum == common::Energy) {
            continue;
        }
        for (const auto& [derivative_key, spectrum_derivative] : flattenedDerivativeSpectra_) {
            const auto& [quantity_enum_derivative, symbol_name] = derivative_key;
            if (quantity_enum_derivative != common::Energy) {
                continue;
            }
            std::pair<common::QuantityEnum, std::pair<common::QuantityEnum, model::symbols::SymbolName>> 
                derivative_product_key = {quantity_enum, derivative_key};

            auto mb_derivative_product = allQuantitiesGetter.getSpectrumDerivativeProduct(
                quantity_enum, quantity_enum_derivative, symbol_name);
            if (mb_derivative_product.has_value()) {
                const auto& spectrum = mb_derivative_product.value();
                flattenedDerivativeProductSpectra_[derivative_product_key] = flatten(spectrum, factories);
            } else {
                flattenedDerivativeProductSpectra_[derivative_product_key] = 
                    transform_one_or_many(
                    std::function([](
                        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& a, 
                        const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& b){
                        return a->element_wise_multiplication(b);
                    }), 
                    spectrum,
                    spectrum_derivative);
            }
        }
    }
}

std::optional<OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>> 
FlattenedSpectra::getFlattenSpectrum(common::QuantityEnum quantity_enum) const {
    if (flattenedSpectra_.contains(quantity_enum)) {
        const auto& spectrum = flattenedSpectra_.at(quantity_enum);
        return copyRef<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>, 
            std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>(spectrum);
    } else {
        return std::nullopt;
    }
}

std::optional<OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>> 
FlattenedSpectra::getFlattenDerivativeSpectrum(common::QuantityEnum quantity_enum, 
        const model::symbols::SymbolName& symbol_name) const {
    if (flattenedDerivativeSpectra_.contains({quantity_enum, symbol_name})) {
        const auto& spectrum = flattenedDerivativeSpectra_.at({quantity_enum, symbol_name});
        return copyRef<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>, 
            std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>(spectrum);
    } else {
        return std::nullopt;
    }
}

std::optional<OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>> 
FlattenedSpectra::getFlattenDerivativeProductSpectrum(common::QuantityEnum quantity_enum, common::QuantityEnum quantity_enum_derivative,
const model::symbols::SymbolName& symbol_name) const {
    std::pair<common::QuantityEnum, std::pair<common::QuantityEnum, model::symbols::SymbolName>> key = {quantity_enum, {quantity_enum_derivative, symbol_name}};
    if (flattenedDerivativeProductSpectra_.contains(key)) {
        const auto& spectrum = flattenedDerivativeProductSpectra_.at(key);
        return copyRef<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>, 
            std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>(spectrum);
    } else {
        return std::nullopt;
    }
}


OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> 
    FlattenedSpectra::getWeights() const {
    return copyRef<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>, 
        std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>(flattenedWeights_);
}


} // namespace eigendecompositor