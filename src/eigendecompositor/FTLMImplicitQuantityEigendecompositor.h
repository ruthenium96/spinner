#ifndef SPINNER_FTLMIMPLICITQUANTITYEIGENDECOMPOSITOR_H
#define SPINNER_FTLMIMPLICITQUANTITYEIGENDECOMPOSITOR_H

#include "ImplicitQuantityEigendecompositor.h"

namespace eigendecompositor {

class FTLMImplicitQuantityEigendecompositor : public ImplicitQuantityEigendecompositor {
  public:
    FTLMImplicitQuantityEigendecompositor(        
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        quantum::linear_algebra::FactoriesList factories_list,
        common::QuantityEnum quantity_implicit_enum,
        uint32_t max_ntz_proj);
    
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        size_t number_of_block,
        const space::Subspace& subspace) override;

    std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
        getSubspectrum(common::QuantityEnum, size_t number_of_block) const override;
    void initialize(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        uint32_t number_of_subspaces) override;
    void finalize() override;

  private:
    Spectrum quantity_implicit_spectrum_multiplied_by_squared_back_projection_;
    common::QuantityEnum quantity_implicit_enum_;
};

} // namespace eigendecompositor

#endif // SPINNER_FTLMIMPLICITQUANTITYEIGENDECOMPOSITOR_H