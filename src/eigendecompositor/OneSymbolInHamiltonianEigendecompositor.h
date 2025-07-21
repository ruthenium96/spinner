#ifndef SPINNER_ONESYMBOLINHAMILTONIANEIGENDECOMPOSITOR_H
#define SPINNER_ONESYMBOLINHAMILTONIANEIGENDECOMPOSITOR_H

#include "AbstractEigendecompositor.h"
#include <functional>

namespace eigendecompositor {

class OneSymbolInHamiltonianEigendecompositor: public AbstractEigendecompositor {
  public:
    OneSymbolInHamiltonianEigendecompositor(
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        std::function<double()> currentValueGetter);
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
    std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
    getSubspectrumDerivativeProduct(common::QuantityEnum, common::QuantityEnum, const model::symbols::SymbolName&, size_t number_of_block) const override;

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
    std::vector<OneOrMany<Subspectrum>> current_energy_spectrum_;
    std::optional<std::vector<OneOrMany<Subspectrum>>> current_energy_derivative_spectrum_;
    std::vector<
        std::optional<OneOrMany<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>>
        eigenvectors_;
#ifndef NDEBUG
    std::vector<OneOrMany<Submatrix>> current_energy_matrix_;
#endif
    double initial_value_of_symbol_;
    bool first_iteration_has_been_done_ = false;
    std::function<double()> currentValueGetter_;
};

}  // namespace eigendecompositor

#endif  //SPINNER_ONESYMBOLINHAMILTONIANEIGENDECOMPOSITOR_H
