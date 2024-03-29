#ifndef SPINNER_IMPLICITSSQUAREEIGENDECOMPOSITOR_H
#define SPINNER_IMPLICITSSQUAREEIGENDECOMPOSITOR_H

#include "AbstractEigendecompositor.h"

namespace eigendecompositor {

// This class implicitly calculate S^2 values using block_properties data (total_mult).
class ImplicitSSquareEigendecompositor: public AbstractEigendecompositor {
  public:
    ImplicitSSquareEigendecompositor(
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        quantum::linear_algebra::FactoriesList factories_list);
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        size_t number_of_block,
        const space::Subspace& subspace) override;

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
    void initialize(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        uint32_t number_of_subspaces) override;
    void finalize() override;

  private:
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor_;
    Spectrum s_square_implicit_spectrum_;
    quantum::linear_algebra::FactoriesList factories_list_;
    bool first_iteration_has_been_done_ = false;
};

}  // namespace eigendecompositor

#endif  //SPINNER_IMPLICITSSQUAREEIGENDECOMPOSITOR_H
