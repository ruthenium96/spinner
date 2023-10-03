#ifndef SPINNER_EXACTEIGENDECOMPOSITOR_H
#define SPINNER_EXACTEIGENDECOMPOSITOR_H

#include <map>
#include <optional>

#include "AbstractEigendecompositor.h"
namespace eigendecompositor {
class ExactEigendecompositor: public AbstractEigendecompositor {
  public:
    void BuildSpectra(
        const std::map<common::QuantityEnum, std::shared_ptr<model::operators::Operator>>&
            operators_,
        const std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<model::operators::Operator>>& derivatives_operators_,
        const space::Space& space,
        const lexicographic::IndexConverter& converter,
        quantum::linear_algebra::FactoriesList data_structure_factories) override;

    const Spectrum& getSpectrum(common::QuantityEnum) const override;
    const Matrix& getMatrix(common::QuantityEnum) const override;
    const Spectrum&
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    const Matrix&
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    void initializeSSquared() override;
    void initializeGSzSquared() override;
    void initializeDerivative(common::QuantityEnum quantity_enum, model::symbols::SymbolName symbol)
        override;

  protected:
    std::map<common::QuantityEnum, common::Quantity> quantities_map_;
    std::map<std::pair<common::QuantityEnum, model::symbols::SymbolName>, common::Quantity>
        derivatives_map_;
    void BuildSpectraWithoutMatrices(
        size_t number_of_blocks,
        const std::map<common::QuantityEnum, std::shared_ptr<model::operators::Operator>>&
            operators_,
        const std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<model::operators::Operator>>& derivatives_operators_,
        const space::Space& space,
        const lexicographic::IndexConverter& converter,
        quantum::linear_algebra::FactoriesList data_structure_factories);
};

}  // namespace eigendecompositor

#endif  //SPINNER_EXACTEIGENDECOMPOSITOR_H
