#include "ExactEigendecompositor.h"

#include <utility>

namespace eigendecompositor {

ExactEigendecompositor::ExactEigendecompositor(
    lexicographic::IndexConverter converter,
    quantum::linear_algebra::FactoriesList factories_list) :
    converter_(std::move(converter)),
    factories_list_(std::move(factories_list)) {}

const Matrix& ExactEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    return quantities_map_.at(quantity_enum).matrix_;
}

const Spectrum& ExactEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    return quantities_map_.at(quantity_enum).spectrum_;
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
    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_,
    const space::Space& space) {
    size_t number_of_blocks = space.getBlocks().size();

    for (const auto& [quanity_enum, _] : operators_) {
        if (!quantities_map_.contains(quanity_enum)) {
            quantities_map_[quanity_enum] = common::Quantity();
        }
    }

    for (const auto& [pair, _] : derivatives_operators_) {
        if (!derivatives_map_.contains(pair)) {
            derivatives_map_[pair] = common::Quantity();
        }
    }

    for (auto& [_, quantity_] : quantities_map_) {
        quantity_.spectrum_.blocks.clear();
        quantity_.spectrum_.blocks.reserve(number_of_blocks);
    }
    for (auto& [_, derivative] : derivatives_map_) {
        derivative.spectrum_.blocks.clear();
        derivative.spectrum_.blocks.reserve(number_of_blocks);
    }
    BuildSpectraWithoutMatrices(number_of_blocks, operators_, derivatives_operators_, space);
}

void ExactEigendecompositor::BuildSpectraWithoutMatrices(
    size_t number_of_blocks,
    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators_,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators_,
    const space::Space& space) {
    for (size_t block = 0; block < number_of_blocks; ++block) {
        std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>
            unitary_transformation_matrix;

        {
            auto& energy = quantities_map_[common::Energy];
            auto hamiltonian_submatrix = Submatrix(
                space.getBlocks()[block],
                *operators_.at(common::Energy),
                converter_,
                factories_list_);
            auto pair = Subspectrum::energy(hamiltonian_submatrix);
            energy.spectrum_.blocks.emplace_back(std::move(pair.first));
            unitary_transformation_matrix = std::move(pair.second);
            energy.matrix_.blocks.emplace_back(std::move(hamiltonian_submatrix));
        }

        for (auto& [quantity_enum, quantity] : quantities_map_) {
            if (quantity_enum == common::Energy) {
                continue;
            }
            auto non_hamiltonian_submatrix = Submatrix(
                space.getBlocks()[block],
                *operators_.at(quantity_enum),
                converter_,
                factories_list_);
            quantity.spectrum_.blocks.emplace_back(
                Subspectrum::non_energy(non_hamiltonian_submatrix, unitary_transformation_matrix));
            quantity.matrix_.blocks.emplace_back(std::move(non_hamiltonian_submatrix));
        }

        for (auto& [pair, derivative] : derivatives_map_) {
            auto derivative_submatrix = Submatrix(
                space.getBlocks()[block],
                *derivatives_operators_.at(pair),
                converter_,
                factories_list_);
            derivative.spectrum_.blocks.emplace_back(
                Subspectrum::non_energy(derivative_submatrix, unitary_transformation_matrix));
            derivative.matrix_.blocks.emplace_back(std::move(derivative_submatrix));
        }
    }
}

}  // namespace eigendecompositor
