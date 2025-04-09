#ifndef SPINNER_ALLQUANTITIESGETTER_H
#define SPINNER_ALLQUANTITIESGETTER_H

#include <optional>
#include <variant>
#include "src/common/Quantity.h"
#include "src/model/symbols/SymbolName.h"

template <typename T>
using OneOrMany = std::variant<T, std::vector<T>>;


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
      virtual std::optional<SpectrumRef>
      getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
      virtual std::optional<MatrixRef>
      getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
};
}

#endif  //SPINNER_ALLQUANTITIESGETTER_H