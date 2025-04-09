#ifndef SPINNER_FTLMEIGENDECOMPOSITOR_H
#define SPINNER_FTLMEIGENDECOMPOSITOR_H

#include <memory>
#include "src/common/Quantity.h"
#include "src/eigendecompositor/ExactEigendecompositor.h"
#include "src/entities/data_structures/AbstractDenseVector.h"
#include "src/entities/spectrum/Spectrum.h"

namespace eigendecompositor {

class FTLMEigendecompositor : public ExactEigendecompositor {
    public:
    FTLMEigendecompositor(
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
        quantum::linear_algebra::FactoriesList factories_list,
        size_t krylov_subspace_size,
        size_t exact_decomposition_threshold,
        size_t number_of_seeds
    );

    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
    BuildSubspectra(
        size_t number_of_block, const space::Subspace& subspace) override;

    std::optional<OneOrMany<std::reference_wrapper<const Subspectrum>>>
        getSubspectrum(common::QuantityEnum, size_t number_of_block) const override;
    std::optional<OneOrMany<std::reference_wrapper<const Submatrix>>>
        getSubmatrix(common::QuantityEnum, size_t number_of_block) const override;
    std::optional<SpectrumRef>
        getSpectrumDerivative(common::QuantityEnum, const model::symbols::SymbolName&) const override;
    std::optional<MatrixRef>
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
    std::vector<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>> seed_vectors_;
    size_t krylov_subspace_size_;
    size_t exact_decomposition_threshold_;
    size_t number_of_seeds_;

    common::Quantity energy_;
    std::shared_ptr<const model::operators::Operator> energy_operator_;
    bool do_we_need_eigenvectors_;
    std::vector<Subspectrum> squared_back_projection_spectrum_;
};

} // namespace eigendecompositor

#endif // SPINNER_FTLMEIGENDECOMPOSITOR_H