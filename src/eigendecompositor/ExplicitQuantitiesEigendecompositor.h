#ifndef SPINNER_EXPLICITQUANTITIESEIGENDECOMPOSITOR_H
#define SPINNER_EXPLICITQUANTITIESEIGENDECOMPOSITOR_H

#include "AbstractEigendecompositor.h"

namespace eigendecompositor {

class ExplicitQuantitiesEigendecompositor: public AbstractEigendecompositor {
  public:
    ExplicitQuantitiesEigendecompositor(
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        quantum::linear_algebra::FactoriesList factories_list);
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

  protected:
    void initialize(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        uint32_t number_of_subspaces) override;
    std::optional<OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
    BuildSubspectra(size_t number_of_block, const space::Subspace& subspace) override;
    void finalize() override;

  private:
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor_;
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter_;
    quantum::linear_algebra::FactoriesList factories_list_;
    std::map<common::QuantityEnum, std::vector<OneOrMany<Subspectrum>>> quantities_spectra_map_;
    std::map<common::QuantityEnum, std::vector<Submatrix>> quantities_matrix_map_;
    std::map<std::pair<common::QuantityEnum, model::symbols::SymbolName>, std::vector<OneOrMany<Subspectrum>>> derivatives_spectra_map_;
    std::map<std::pair<common::QuantityEnum, model::symbols::SymbolName>, std::vector<Submatrix>> derivatives_matrix_map_;
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>
        quantities_operators_map_;
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>
        derivatives_operators_map_;

    static OneOrMany<Subspectrum> non_energy_subspectrum(
        const Submatrix& non_hamiltonian_submatrix,
        const OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>&
            unitary_transformation_matrix);
};

}  // namespace eigendecompositor

#endif  //SPINNER_EXPLICITQUANTITIESEIGENDECOMPOSITOR_H
