#ifndef SPINNER_ABSTRACTEIGENDECOMPOSITOR_H
#define SPINNER_ABSTRACTEIGENDECOMPOSITOR_H

#include "src/eigendecompositor/AllQuantitiesGetter.h"

namespace eigendecompositor {

class AbstractEigendecompositor : public AllQuantitiesGetter {
  public:
    void BuildSpectra(
        const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators,
        const std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators,
        const space::Space& space);
    bool BuildSpectraWasCalled() const;

    std::optional<OneOrMany<SpectrumRef>> getSpectrum(common::QuantityEnum) const override;
    std::optional<OneOrMany<MatrixRef>> getMatrix(common::QuantityEnum) const override;
    std::optional<OneOrMany<SpectrumRef>>
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    std::optional<OneOrMany<MatrixRef>>
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    size_t getSubspectrumSize(common::QuantityEnum, size_t number_of_block) const;

    virtual ~AbstractEigendecompositor() = default;

    virtual void initialize(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        uint32_t number_of_subspaces) = 0;
    virtual std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(size_t number_of_block, const space::Subspace& subspace) = 0;
    virtual void finalize() = 0;
  private:
    bool buildSpectraWasCalled = false;
    uint32_t number_of_subspaces_;
};

}  // namespace eigendecompositor

#endif  //SPINNER_ABSTRACTEIGENDECOMPOSITOR_H
