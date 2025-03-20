#include "FlattenedSpectra.h"

#include <memory>
#include <optional>

#include "magic_enum.hpp"
#include "src/common/Quantity.h"

namespace eigendecompositor {

void FlattenedSpectra::updateValues(const AllQuantitiesGetter& allQuantitiesGetter,
    const quantum::linear_algebra::FactoriesList& factories) {
    flattenedSpectra_.clear();
    degeneracyValues_ = factories.createVector();
    for (const auto& quantity_enum : magic_enum::enum_values<common::QuantityEnum>()) {
        auto maybe_spectrum = allQuantitiesGetter.getSpectrum(quantity_enum);
        if (maybe_spectrum.has_value()) {
            flattenedSpectra_[quantity_enum] = factories.createVector();
            for (const auto& subspectrum : maybe_spectrum.value().get().blocks) {
                flattenedSpectra_[quantity_enum]->concatenate_with(subspectrum.raw_data);
                if (quantity_enum == common::Energy) {
                    degeneracyValues_->add_identical_values(
                        subspectrum.raw_data->size(),
                        subspectrum.properties.degeneracy);
                }
            }
        }
    }
    flattenedSpectra_[common::Energy]->subtract_minimum();
}

void FlattenedSpectra::updateDerivativeValues(const AllQuantitiesGetter& allQuantitiesGetter,
    const std::vector<model::symbols::SymbolName>& symbol_names,
    const quantum::linear_algebra::FactoriesList& factories) {
    flattenedDerivativeSpectra_.clear();
    for (const auto& quantity_enum : magic_enum::enum_values<common::QuantityEnum>()) {
        for (const auto& symbol_name : symbol_names) {
            auto maybe_spectrum = allQuantitiesGetter.getSpectrumDerivative(quantity_enum, symbol_name);
            if (maybe_spectrum.has_value()) {
                flattenedDerivativeSpectra_[{quantity_enum, symbol_name}] = factories.createVector();
                for (const auto& subspectrum : maybe_spectrum.value().get().blocks) {
                    flattenedDerivativeSpectra_[{quantity_enum, symbol_name}]->concatenate_with(subspectrum.raw_data);
                }
            }    
        }
    }   
}

std::optional<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> 
    FlattenedSpectra::getFlattenSpectrum(common::QuantityEnum quantity_enum) const {
    if (flattenedSpectra_.contains(quantity_enum)) {
        return flattenedSpectra_.at(quantity_enum);
    } else {
        return std::nullopt;
    }
}

std::optional<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> 
    FlattenedSpectra::getFlattenDerivativeSpectrum(common::QuantityEnum quantity_enum, 
        const model::symbols::SymbolName& symbol_name) const {
    if (flattenedDerivativeSpectra_.contains({quantity_enum, symbol_name})) {
        return flattenedDerivativeSpectra_.at({quantity_enum, symbol_name});
    } else {
        return std::nullopt;
    }
}


const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& 
    FlattenedSpectra::getDegeneracyValues() const {
    return degeneracyValues_;
}


} // namespace eigendecompositor