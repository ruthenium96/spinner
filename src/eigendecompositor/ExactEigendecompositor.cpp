#include "ExactEigendecompositor.h"

#include <cassert>
#include <utility>

namespace eigendecompositor {

ExactEigendecompositor::ExactEigendecompositor(
    lexicographic::IndexConverter converter,
    quantum::linear_algebra::FactoriesList factories_list) :
    converter_(std::move(converter)),
    factories_list_(std::move(factories_list)) {}

std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
ExactEigendecompositor::BuildSubspectra(
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_to_calculate,
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_to_calculate,
    size_t number_of_block,
    const space::Subspace& subspace) {
    std::optional<std::shared_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>>
        mb_unitary_transformation_matrix;

    auto hamiltonian_submatrix = Submatrix(
        subspace,
        *operators_to_calculate.at(common::Energy),
        converter_,
        factories_list_);

    if (operators_to_calculate.size() == 1 && derivatives_operators_to_calculate.empty()) {
        // if we need to explicitly calculate _only_ energy, we do not need eigenvectors:
        auto energy_spectrum = energy_subspectrum_eigenvalues_only(hamiltonian_submatrix);
        energy_.spectrum_.blocks.emplace_back(std::move(energy_spectrum));
    } else {
        auto pair = energy_subspectrum_with_eigenvectors(hamiltonian_submatrix);
        mb_unitary_transformation_matrix = std::move(pair.second);
        energy_.spectrum_.blocks.emplace_back(std::move(pair.first));
    }
    energy_.matrix_.blocks.emplace_back(std::move(hamiltonian_submatrix));
    assert(energy_.spectrum_.blocks.size() == number_of_block + 1);
    assert(energy_.matrix_.blocks.size() == number_of_block + 1);

    // delete energy operator, because we just have calculated energy
    std::erase_if(operators_to_calculate, [](auto p) { return p.first == common::Energy; });

    return mb_unitary_transformation_matrix;
}

void ExactEigendecompositor::initialize() {
    energy_.matrix_.blocks.clear();
    energy_.spectrum_.blocks.clear();
}

void ExactEigendecompositor::finalize() {}

std::optional<std::reference_wrapper<const Matrix>>
ExactEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::Energy) {
        return energy_.matrix_;
    }
    return std::nullopt;
}

std::optional<std::reference_wrapper<const Spectrum>>
ExactEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::Energy) {
        return energy_.spectrum_;
    }
    return std::nullopt;
}

std::optional<std::reference_wrapper<const Spectrum>> ExactEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return std::nullopt;
}

std::optional<std::reference_wrapper<const Matrix>> ExactEigendecompositor::getMatrixDerivative(
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
