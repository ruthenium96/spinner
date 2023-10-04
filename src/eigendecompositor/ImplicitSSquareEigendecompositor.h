#ifndef SPINNER_IMPLICITSSQUAREEIGENDECOMPOSITOR_H
#define SPINNER_IMPLICITSSQUAREEIGENDECOMPOSITOR_H

#include "AbstractEigendecompositor.h"

namespace eigendecompositor {

class ImplicitSSquareEigendecompositor: public AbstractEigendecompositor {
  public:
    ImplicitSSquareEigendecompositor(
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        quantum::linear_algebra::FactoriesList factories_list);
    void BuildSpectra(
        const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_,
        const std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_,
        const space::Space& space) override;

    std::optional<std::reference_wrapper<const Spectrum>>
    getSpectrum(common::QuantityEnum anEnum) const override;
    std::optional<std::reference_wrapper<const Matrix>>
    getMatrix(common::QuantityEnum anEnum) const override;
    std::optional<std::reference_wrapper<const Spectrum>> getSpectrumDerivative(
        common::QuantityEnum anEnum,
        const model::symbols::SymbolName& name) const override;
    std::optional<std::reference_wrapper<const Matrix>> getMatrixDerivative(
        common::QuantityEnum anEnum,
        const model::symbols::SymbolName& name) const override;

  private:
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor_;
    common::Quantity s_square_implicit_;
    quantum::linear_algebra::FactoriesList factories_list_;
};

}  // namespace eigendecompositor

#endif  //SPINNER_IMPLICITSSQUAREEIGENDECOMPOSITOR_H
