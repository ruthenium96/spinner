#ifndef SPINNER_EXPLICITQUANTITIESEIGENDECOMPOSITOR_H
#define SPINNER_EXPLICITQUANTITIESEIGENDECOMPOSITOR_H

#include "AbstractEigendecompositor.h"

namespace eigendecompositor {

class ExplicitQuantitiesEigendecompositor: public AbstractEigendecompositor {
  public:
    ExplicitQuantitiesEigendecompositor(
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        lexicographic::IndexConverter converter,
        quantum::linear_algebra::FactoriesList factories_list);
    std::optional<std::reference_wrapper<const Spectrum>>
    getSpectrum(common::QuantityEnum quantity_enum) const override;
    std::optional<std::reference_wrapper<const Matrix>>
    getMatrix(common::QuantityEnum quantity_enum) const override;
    std::optional<std::reference_wrapper<const Spectrum>> getSpectrumDerivative(
        common::QuantityEnum quantity_enum,
        const model::symbols::SymbolName& symbol_name) const override;
    std::optional<std::reference_wrapper<const Matrix>> getMatrixDerivative(
        common::QuantityEnum quantity_enum,
        const model::symbols::SymbolName& symbol_name) const override;

  protected:
    void initialize() override;
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        size_t number_of_block,
        const space::Subspace& subspace) override;
    void finalize() override;

  private:
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor_;
    lexicographic::IndexConverter converter_;
    quantum::linear_algebra::FactoriesList factories_list_;
    std::map<common::QuantityEnum, common::Quantity> quantities_map_;
    std::map<std::pair<common::QuantityEnum, model::symbols::SymbolName>, common::Quantity>
        derivatives_map_;

    static Subspectrum non_energy_subspectrum(
        const Submatrix& non_hamiltonian_submatrix,
        const std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>&
            unitary_transformation_matrix);
};

}  // namespace eigendecompositor

#endif  //SPINNER_EXPLICITQUANTITIESEIGENDECOMPOSITOR_H
