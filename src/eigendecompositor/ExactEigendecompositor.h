#ifndef SPINNER_EXACTEIGENDECOMPOSITOR_H
#define SPINNER_EXACTEIGENDECOMPOSITOR_H

#include <map>
#include <optional>

#include "AbstractEigendecompositor.h"
namespace eigendecompositor {
class ExactEigendecompositor: public AbstractEigendecompositor {
  public:
    ExactEigendecompositor(
        lexicographic::IndexConverter converter,
        quantum::linear_algebra::FactoriesList factories_list);
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_,
        size_t number_of_block,
        const space::Subspace& subspace) override;

    std::optional<std::reference_wrapper<const Spectrum>>
        getSpectrum(common::QuantityEnum) const override;
    std::optional<std::reference_wrapper<const Matrix>>
        getMatrix(common::QuantityEnum) const override;
    std::optional<std::reference_wrapper<const Spectrum>>
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    std::optional<std::reference_wrapper<const Matrix>>
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    void initialize() override;
    void finalize() override;

  private:
    lexicographic::IndexConverter converter_;
    quantum::linear_algebra::FactoriesList factories_list_;
    common::Quantity energy_;

    static Subspectrum energy_subspectrum_eigenvalues_only(const Submatrix& hamiltonian_submatrix);
    static std::
        pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        energy_subspectrum_with_eigenvectors(const Submatrix& hamiltonian_submatrix);
};

}  // namespace eigendecompositor

#endif  //SPINNER_EXACTEIGENDECOMPOSITOR_H
