#ifndef SPINNER_FLATTENEDSPECTRA_H
#define SPINNER_FLATTENEDSPECTRA_H

#include <functional>
#include <map>
#include <memory>

#include "src/common/Quantity.h"
#include "src/eigendecompositor/AllQuantitiesGetter.h"
#include "src/entities/data_structures/AbstractDenseVector.h"

namespace eigendecompositor {
class FlattenedSpectra {
  public:
    FlattenedSpectra() = default;
    void updateValues(const AllQuantitiesGetter& allQuantitiesGetter,
        const quantum::linear_algebra::FactoriesList& factories);
    std::optional<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>> 
        getFlattenSpectrum(common::QuantityEnum quantity_enum) const;
    const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>& getDegeneracyValues() const;

  private:
    std::map<common::QuantityEnum, 
        std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> flattenedSpectra_;
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> degeneracyValues_;
};
} // namespace eigendecompositor

#endif // SPINNER_FLATTENEDSPECTRA_H