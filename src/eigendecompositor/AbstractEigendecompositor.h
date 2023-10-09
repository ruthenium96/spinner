#ifndef SPINNER_ABSTRACTEIGENDECOMPOSITOR_H
#define SPINNER_ABSTRACTEIGENDECOMPOSITOR_H

#include <optional>

#include "src/model/Model.h"

namespace eigendecompositor {
class AbstractEigendecompositor {
  public:
    void BuildSpectra(
        const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators,
        const std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators,
        const space::Space& space);

    virtual std::optional<std::reference_wrapper<const Spectrum>>
        getSpectrum(common::QuantityEnum) const = 0;
    virtual std::optional<std::reference_wrapper<const Matrix>>
        getMatrix(common::QuantityEnum) const = 0;
    virtual std::optional<std::reference_wrapper<const Spectrum>>
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;
    virtual std::optional<std::reference_wrapper<const Matrix>>
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const = 0;

    virtual ~AbstractEigendecompositor() = default;

    virtual void initialize() = 0;
    virtual std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        size_t number_of_block,
        const space::Subspace& subspace) = 0;
    virtual void finalize() = 0;
};

}  // namespace eigendecompositor

#endif  //SPINNER_ABSTRACTEIGENDECOMPOSITOR_H
