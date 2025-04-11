#ifndef SPINNER_ALLQUANTITIESGETTER_H
#define SPINNER_ALLQUANTITIESGETTER_H

#include <optional>

#include "src/common/OneOrMany.h"
#include "src/common/Quantity.h"
#include "src/model/symbols/SymbolName.h"

namespace eigendecompositor {
class AllQuantitiesGetter {
    public:
      // What means Spectrum's optional?
      // The getter *may* return std::nullopt if corresponding Operator has been passed into initialize.
      // What means Matrix's optional?
      // The getter returns some value if underlying wrapper actually store Matrix.
      virtual std::optional<OneOrMany<SpectrumRef>>
          getSpectrum(common::QuantityEnum) const = 0;
      virtual std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
          getSubspectrum(common::QuantityEnum, size_t number_of_block) const = 0;
      virtual std::optional<OneOrMany<MatrixRef>>
          getMatrix(common::QuantityEnum) const = 0;
      virtual std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
          getSubmatrix(common::QuantityEnum, size_t number_of_block) const = 0;
      virtual std::optional<OneOrMany<SpectrumRef>>
      getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
      virtual std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
      getSubspectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&, size_t number_of_block) const = 0;
      virtual std::optional<OneOrMany<MatrixRef>>
      getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
      virtual std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
      getSubmatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&, size_t number_of_block) const = 0;
};
}

#endif  //SPINNER_ALLQUANTITIESGETTER_H