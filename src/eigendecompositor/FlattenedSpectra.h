#ifndef SPINNER_FLATTENEDSPECTRA_H
#define SPINNER_FLATTENEDSPECTRA_H

#include <functional>
#include <map>
#include <memory>

#include "src/common/Quantity.h"
#include "src/eigendecompositor/AllQuantitiesGetter.h"
#include "src/linalg_structures/AbstractDenseVector.h"
#include "src/model/symbols/SymbolName.h"

namespace spinner::eigendecompositor {
class FlattenedSpectra {
  public:
    FlattenedSpectra() = default;
    void updateValues(const AllQuantitiesGetter& allQuantitiesGetter,
        const linalg_structures::FactoriesList& factories);
    void updateDerivativeValues(const AllQuantitiesGetter& allQuantitiesGetter,
        const std::vector<model::symbols::SymbolName>& symbol_names,
        const linalg_structures::FactoriesList& factories);
  
    std::optional<OneOrMany<std::reference_wrapper<const std::unique_ptr<linalg_structures::AbstractDenseVector>>>> 
        getFlattenSpectrum(common::QuantityEnum quantity_enum) const;
    std::optional<OneOrMany<std::reference_wrapper<const std::unique_ptr<linalg_structures::AbstractDenseVector>>>> 
        getFlattenDerivativeSpectrum(common::QuantityEnum quantity_enum, 
          const model::symbols::SymbolName& symbol_name) const;
    std::optional<OneOrMany<std::reference_wrapper<const std::unique_ptr<linalg_structures::AbstractDenseVector>>>> 
        getFlattenDerivativeProductSpectrum(common::QuantityEnum quantity_enum, common::QuantityEnum quantity_enum_derivative,
        const model::symbols::SymbolName& symbol_name) const;
  
    OneOrMany<std::reference_wrapper<const std::unique_ptr<linalg_structures::AbstractDenseVector>>> getWeights() const;

  private:
    std::map<common::QuantityEnum, 
        OneOrMany<std::unique_ptr<linalg_structures::AbstractDenseVector>>> flattenedSpectra_;
    std::map<std::pair<common::QuantityEnum, model::symbols::SymbolName>, 
        OneOrMany<std::unique_ptr<linalg_structures::AbstractDenseVector>>> flattenedDerivativeSpectra_;
    std::map<std::pair<common::QuantityEnum, std::pair<common::QuantityEnum, model::symbols::SymbolName>>, 
        OneOrMany<std::unique_ptr<linalg_structures::AbstractDenseVector>>> flattenedDerivativeProductSpectra_;

    std::unique_ptr<linalg_structures::AbstractDenseVector> degeneracyValues_;
    OneOrMany<std::unique_ptr<linalg_structures::AbstractDenseVector>> flattenedWeights_;
};
} // namespace spinner::eigendecompositor

#endif // SPINNER_FLATTENEDSPECTRA_H