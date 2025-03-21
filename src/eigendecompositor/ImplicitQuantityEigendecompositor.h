#ifndef SPINNER_IMPLICITQUANTITYEIGENDECOMPOSITOR_H
#define SPINNER_IMPLICITQUANTITYEIGENDECOMPOSITOR_H

#include "AbstractEigendecompositor.h"
#include "src/common/Quantity.h"

namespace eigendecompositor {

// This class implicitly calculate S^2 or M^2 values using block_properties data (total_mult or n_proj).
class ImplicitQuantityEigendecompositor: public AbstractEigendecompositor {
  public:
    ImplicitQuantityEigendecompositor(
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        quantum::linear_algebra::FactoriesList factories_list,
        common::QuantityEnum quantity_implicit_enum,
        uint32_t max_ntz_proj);
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        size_t number_of_block,
        const space::Subspace& subspace) override;

    std::optional<SpectrumRef>
    getSpectrum(common::QuantityEnum anEnum) const override;
    std::optional<MatrixRef>
    getMatrix(common::QuantityEnum anEnum) const override;
    std::optional<SpectrumRef> getSpectrumDerivative(
        common::QuantityEnum anEnum,
        const model::symbols::SymbolName& name) const override;
    std::optional<MatrixRef> getMatrixDerivative(
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
    uint32_t max_ntz_proj_;
    Spectrum quantity_implicit_spectrum_;
    common::QuantityEnum quantity_implicit_enum_;
    quantum::linear_algebra::FactoriesList factories_list_;
    bool first_iteration_has_been_done_ = false;

    double calculate_value(const BlockProperties& properties) const;
};

}  // namespace eigendecompositor

#endif  //SPINNER_IMPLICITQUANTITYEIGENDECOMPOSITOR_H
