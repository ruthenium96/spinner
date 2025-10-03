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
    std::optional<OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
    BuildSubspectra(
        size_t number_of_block,
        const space::Subspace& subspace) override;

    std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
    getSubspectrum(common::QuantityEnum, size_t number_of_block) const override;
    std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
    getSubmatrix(common::QuantityEnum, size_t number_of_block) const override;
    std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
    getSubspectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&, size_t number_of_block) const override;
    std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
    getSubmatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&, size_t number_of_block) const override;

    OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>
    getWeightsOfBlockStates(size_t number_of_block) const override;

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
