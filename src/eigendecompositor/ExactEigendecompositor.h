#ifndef SPINNER_EXACTEIGENDECOMPOSITOR_H
#define SPINNER_EXACTEIGENDECOMPOSITOR_H

#include <map>
#include <optional>

#include "AbstractEigendecompositor.h"

namespace eigendecompositor {
class ExactEigendecompositor: public AbstractEigendecompositor {
  public:
    ExactEigendecompositor(
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        quantum::linear_algebra::FactoriesList factories_list);
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        size_t number_of_block, const space::Subspace& subspace) override;

    std::optional<std::reference_wrapper<const Spectrum>>
        getSpectrum(common::QuantityEnum) const override;
    std::optional<std::reference_wrapper<const Matrix>>
        getMatrix(common::QuantityEnum) const override;
    std::optional<std::reference_wrapper<const Spectrum>>
    getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    std::optional<std::reference_wrapper<const Matrix>>
    getMatrixDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    void initialize(
        std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
            operators_to_calculate,
        std::map<
            std::pair<common::QuantityEnum, model::symbols::SymbolName>,
            std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
        uint32_t number_of_subspaces) override;
    void finalize() override;

  private:
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter_;
    quantum::linear_algebra::FactoriesList factories_list_;
    common::Quantity energy_;
    std::shared_ptr<const model::operators::Operator> energy_operator_;
    bool do_we_need_eigenvectors_;

    static Subspectrum energy_subspectrum_eigenvalues_only(const Submatrix& hamiltonian_submatrix);
    static std::
        pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        energy_subspectrum_with_eigenvectors(const Submatrix& hamiltonian_submatrix);
};

}  // namespace eigendecompositor

#endif  //SPINNER_EXACTEIGENDECOMPOSITOR_H
