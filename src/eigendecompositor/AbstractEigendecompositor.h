#ifndef SPINNER_ABSTRACTEIGENDECOMPOSITOR_H
#define SPINNER_ABSTRACTEIGENDECOMPOSITOR_H

#include <optional>

#include "src/model/Model.h"

namespace eigendecompositor {
class AbstractEigendecompositor {
  public:
    virtual void BuildSpectra(
        const std::map<common::QuantityEnum, std::shared_ptr<model::operators::Operator>>&
            operators_,
        const std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<model::operators::Operator>>& derivatives_operators_,
        const space::Space& space,
        const lexicographic::IndexConverter& converter,
        quantum::linear_algebra::FactoriesList data_structure_factories) = 0;

    virtual const Spectrum& getSpectrum(common::QuantityEnum) const = 0;
    virtual const Matrix& getMatrix(common::QuantityEnum) const = 0;
    virtual const Spectrum&
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
    virtual const Matrix&
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;

    virtual ~AbstractEigendecompositor() = default;
};

}  // namespace eigendecompositor

#endif  //SPINNER_ABSTRACTEIGENDECOMPOSITOR_H
