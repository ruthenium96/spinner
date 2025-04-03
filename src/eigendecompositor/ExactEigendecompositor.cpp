#include "ExactEigendecompositor.h"

#include <functional>
#include <utility>
#include <vector>

namespace eigendecompositor {

ExactEigendecompositor::ExactEigendecompositor(
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter,
    quantum::linear_algebra::FactoriesList factories_list) :
    converter_(std::move(converter)),
    factories_list_(std::move(factories_list)) {}

std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
ExactEigendecompositor::BuildSubspectra(
    size_t number_of_block,
    const space::Subspace& subspace) {
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        mb_unitary_transformation_matrix;

    // return_sparse_if_possible is false, because eigendecomposition of dense matrix is faster
    auto hamiltonian_submatrix = Submatrix(subspace, *energy_operator_, converter_, factories_list_, false);

    if (!do_we_need_eigenvectors_) {
        // if we need to explicitly calculate _only_ energy, we do not need eigenvectors:
        auto energy_spectrum = energy_subspectrum_eigenvalues_only(hamiltonian_submatrix);
        energy_.spectrum_.blocks[number_of_block] = std::move(energy_spectrum);
    } else {
        auto pair = energy_subspectrum_with_eigenvectors(hamiltonian_submatrix);
        mb_unitary_transformation_matrix = std::move(pair.second);
        energy_.spectrum_.blocks[number_of_block] = std::move(pair.first);
    }
#ifndef NDEBUG
    energy_.matrix_.blocks[number_of_block] = std::move(hamiltonian_submatrix);
#endif

    return mb_unitary_transformation_matrix;
}

void ExactEigendecompositor::initialize(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    uint32_t number_of_subspaces) {
    energy_.matrix_.blocks.clear();
    energy_.spectrum_.blocks.clear();
    energy_.matrix_.blocks.resize(number_of_subspaces);
    energy_.spectrum_.blocks.resize(number_of_subspaces);

    energy_operator_ = operators_to_calculate.at(common::Energy);
    do_we_need_eigenvectors_ =
        !(operators_to_calculate.size() == 1 && derivatives_operators_to_calculate.empty());
    std::erase_if(operators_to_calculate, [](const auto& p) { return p.first == common::Energy; });
}

void ExactEigendecompositor::finalize() {}

std::optional<MatrixRef>
ExactEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
#ifndef NDEBUG
    if (quantity_enum == common::Energy) {
        return MatrixRef(energy_.matrix_);
    }
#endif
    return std::nullopt;
}

std::optional<SpectrumRef>
ExactEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::Energy) {
        return SpectrumRef(energy_.spectrum_);
    }
    return std::nullopt;
}

std::optional<std::reference_wrapper<const Subspectrum>>
ExactEigendecompositor::getSubspectrum(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    if (quantity_enum == common::Energy) {
        return energy_.spectrum_.blocks[number_of_block];
    }
    return std::nullopt;
}

std::optional<SpectrumRef> ExactEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return std::nullopt;
}

std::optional<MatrixRef> ExactEigendecompositor::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return std::nullopt;
}

Subspectrum ExactEigendecompositor::energy_subspectrum_eigenvalues_only(
    const Submatrix& hamiltonian_submatrix) {
    auto eigenvalues = hamiltonian_submatrix.raw_data->diagonalizeValues();

    auto energy_subspectrum = Subspectrum(std::move(eigenvalues), hamiltonian_submatrix.properties);

    return std::move(energy_subspectrum);
}

std::pair<Subspectrum, std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
ExactEigendecompositor::energy_subspectrum_with_eigenvectors(
    const Submatrix& hamiltonian_submatrix) {
    auto eigencouple = hamiltonian_submatrix.raw_data->diagonalizeValuesVectors();

    auto energy_subspectrum =
        Subspectrum(std::move(eigencouple.eigenvalues), hamiltonian_submatrix.properties);

    return {std::move(energy_subspectrum), std::move(eigencouple.eigenvectors)};
}

}  // namespace eigendecompositor
