#include "ExactEigendecompositor.h"

#include <cassert>

namespace eigendecompositor {

const Matrix& ExactEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return energy.matrix_;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return s_squared->matrix_;
    } else if (quantity_enum == common::QuantityEnum::gSz_total_squared) {
        return g_sz_squared->matrix_;
    }
    throw std::invalid_argument("There is no such matrix");
    assert(0);
}

const Spectrum& ExactEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return energy.spectrum_;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return s_squared->spectrum_;
    } else if (quantity_enum == common::QuantityEnum::gSz_total_squared) {
        return g_sz_squared->spectrum_;
    }
    throw std::invalid_argument("There is no such spectrum");
    assert(0);
}

const Spectrum& ExactEigendecompositor::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return derivatives_map_.at({quantity_enum, symbol}).spectrum_;
}

const Matrix& ExactEigendecompositor::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return derivatives_map_.at({quantity_enum, symbol}).matrix_;
}

void ExactEigendecompositor::BuildSpectra(
    std::shared_ptr<const model::operators::Operator> energy_operator,
    std::optional<std::shared_ptr<const model::operators::Operator>> s_squared_operator,
    std::optional<std::shared_ptr<const model::operators::Operator>> g_sz_squared_operator,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<model::operators::Operator>>& derivatives_operators_,
    const space::Space& space,
    const lexicographic::IndexConverter& converter,
    quantum::linear_algebra::FactoriesList data_structure_factories) {
    size_t number_of_blocks = space.getBlocks().size();

    energy.spectrum_.blocks.clear();
    energy.spectrum_.blocks.reserve(number_of_blocks);
    if (s_squared.has_value()) {
        s_squared->spectrum_.blocks.clear();
        s_squared->spectrum_.blocks.reserve(number_of_blocks);
    }
    if (g_sz_squared.has_value()) {
        g_sz_squared->spectrum_.blocks.clear();
        g_sz_squared->spectrum_.blocks.reserve(number_of_blocks);
    }
    for (auto& [_, derivative] : derivatives_map_) {
        derivative.spectrum_.blocks.clear();
        derivative.spectrum_.blocks.reserve(number_of_blocks);
    }
    BuildSpectraWithoutMatrices(
        number_of_blocks,
        energy_operator,
        s_squared_operator,
        g_sz_squared_operator,
        derivatives_operators_,
        space,
        converter,
        data_structure_factories);
}

void ExactEigendecompositor::BuildSpectraWithoutMatrices(
    size_t number_of_blocks,
    std::shared_ptr<const model::operators::Operator> energy_operator,
    std::optional<std::shared_ptr<const model::operators::Operator>> s_squared_operator,
    std::optional<std::shared_ptr<const model::operators::Operator>> g_sz_squared_operator,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<model::operators::Operator>>& derivatives_operators_,
    const space::Space& space,
    const lexicographic::IndexConverter& converter,
    quantum::linear_algebra::FactoriesList data_structure_factories) {
    for (size_t block = 0; block < number_of_blocks; ++block) {
        std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>
            unitary_transformation_matrix;

        {
            auto hamiltonian_submatrix = Submatrix(
                space.getBlocks()[block],
                *energy_operator,
                converter,
                data_structure_factories);
            auto pair = Subspectrum::energy(hamiltonian_submatrix);
            energy.spectrum_.blocks.emplace_back(std::move(pair.first));
            unitary_transformation_matrix = std::move(pair.second);
            energy.matrix_.blocks.emplace_back(std::move(hamiltonian_submatrix));
        }

        if (s_squared.has_value()) {
            auto non_hamiltonian_submatrix = Submatrix(
                space.getBlocks()[block],
                *s_squared_operator.value(),
                converter,
                data_structure_factories);
            s_squared->spectrum_.blocks.emplace_back(
                Subspectrum::non_energy(non_hamiltonian_submatrix, unitary_transformation_matrix));
            s_squared->matrix_.blocks.emplace_back(std::move(non_hamiltonian_submatrix));
        }

        if (g_sz_squared.has_value()) {
            auto non_hamiltonian_submatrix = Submatrix(
                space.getBlocks()[block],
                *g_sz_squared_operator.value(),
                converter,
                data_structure_factories);
            g_sz_squared->spectrum_.blocks.emplace_back(
                Subspectrum::non_energy(non_hamiltonian_submatrix, unitary_transformation_matrix));
            g_sz_squared->matrix_.blocks.emplace_back(std::move(non_hamiltonian_submatrix));
        }

        for (auto& [pair, derivative] : derivatives_map_) {
            auto derivative_submatrix = Submatrix(
                space.getBlocks()[block],
                *derivatives_operators_.at(pair),
                converter,
                data_structure_factories);
            derivative.spectrum_.blocks.emplace_back(
                Subspectrum::non_energy(derivative_submatrix, unitary_transformation_matrix));
            derivative.matrix_.blocks.emplace_back(std::move(derivative_submatrix));
        }
    }
}

void ExactEigendecompositor::initializeSSquared() {
    s_squared = common::Quantity();
}

void ExactEigendecompositor::initializeGSzSquared() {
    g_sz_squared = common::Quantity();
}

void ExactEigendecompositor::initializeDerivative(
    common::QuantityEnum quantity_enum,
    model::symbols::SymbolName symbol) {
    derivatives_map_[{quantity_enum, symbol}] = common::Quantity();
}

}  // namespace eigendecompositor
