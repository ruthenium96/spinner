#ifndef SPINNER_ALLQUANTITIESGETTER_H
#define SPINNER_ALLQUANTITIESGETTER_H

#include <optional>
#include "src/common/Quantity.h"
#include "src/model/symbols/SymbolName.h"

namespace eigendecompositor {
class AllQuantitiesGetter {
    public:
      // What means Spectrum's optional?
      // The getter *may* return std::nullopt if corresponding Operator has been passed into initialize.
      // What means Matrix's optional?
      // The getter returns some value if underlying wrapper actually store Matrix.
      virtual std::optional<std::reference_wrapper<const Spectrum>>
          getSpectrum(common::QuantityEnum) const = 0;
      virtual std::optional<std::reference_wrapper<const Matrix>>
          getMatrix(common::QuantityEnum) const = 0;
      virtual std::optional<std::reference_wrapper<const Spectrum>>
      getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
      virtual std::optional<std::reference_wrapper<const Matrix>>
      getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
};
}

#endif  //SPINNER_ALLQUANTITIESGETTER_H