#ifndef SPINNER_ONESYMBOLINHAMILTONIANEIGENDECOMPOSITOR_H
#define SPINNER_ONESYMBOLINHAMILTONIANEIGENDECOMPOSITOR_H

#include "AbstractEigendecompositor.h"

namespace eigendecompositor {

class OneSymbolInHamiltonianEigendecompositor: public AbstractEigendecompositor {
  public:
    OneSymbolInHamiltonianEigendecompositor(
        std::unique_ptr<AbstractEigendecompositor> eigendecompositor,
        std::function<double()> currentValueGetter);
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        size_t number_of_block,
        const space::Subspace& subspace) override;

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
    void initialize() override;
    void finalize() override;

  private:
    std::unique_ptr<AbstractEigendecompositor> eigendecompositor_;
    common::Quantity current_energy_;
    std::optional<common::Quantity> current_energy_derivative_;
    std::vector<
        std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>>
        eigenvectors_;
    double initial_value_of_symbol_;
    bool first_iteration_has_been_done_ = false;
    std::function<double()> currentValueGetter_;
};

}  // namespace eigendecompositor

#endif  //SPINNER_ONESYMBOLINHAMILTONIANEIGENDECOMPOSITOR_H
