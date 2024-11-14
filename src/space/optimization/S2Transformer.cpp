#include "S2Transformer.h"

#include "src/common/Logger.h"
#include <utility>

namespace space::optimization {

S2Transformer::S2Transformer(
    std::shared_ptr<const lexicographic::IndexConverter> converter,
    quantum::linear_algebra::FactoriesList factories,
    std::shared_ptr<const spin_algebra::SSquaredConverter> ssquared_converter,
    bool ito) :
    ito_(ito),
    converter_(std::move(converter)),
    factories_(std::move(factories)),
    ssquared_converter_(std::move(ssquared_converter)) {}

space::Space S2Transformer::apply(Space&& space) const {
    std::vector<Subspace> vector_result;
    vector_result.reserve(space.getBlocks().size());

    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace = space.getBlocks()[i];

        // It is total_proj + max_total_spin = ntz_value
        // We need 2 * |total_proj| + 1 = 2 * |ntz_value - max_total_spin| + 1 =
        // |2 * ntz_value - 2 * max_total_spin - 1 + 1| + 1 =
        // |2 * ntz_value - max_total_proj + 1| + 1
        auto ntz_value = subspace.properties.n_proj.value();
        int a = 2 * (int)ntz_value + 1 - (int)converter_->get_max_ntz_proj();
        spin_algebra::Multiplicity current_mult = std::abs(a) + 1;

        spin_algebra::SSquaredState::Properties subspace_properties;
        subspace_properties.multiplicity = current_mult;
        subspace_properties.representations = subspace.properties.representation;

        if (ito_) {
            auto mb_indexes_of_block = ssquared_converter_->indexes_with_property(subspace_properties);
            if (!mb_indexes_of_block.has_value()) {
                continue;
            }
            const auto& indexes_of_block = mb_indexes_of_block.value();

            subspace.ssquared_indexes = indexes_of_block;
        } else {
            auto mb_ssquared_states = ssquared_converter_->block_with_property(subspace_properties);
            if (!mb_ssquared_states.has_value()) {
                continue;
            }
            const auto& s_squared_states = mb_ssquared_states.value().get();

            subspace.dense_semiunitary_matrix =
                constructTransformationMatrix(s_squared_states, subspace);
        }

        subspace.properties.total_mult = current_mult;
        subspace.properties.degeneracy *= current_mult;
        // TODO: for cases with TzSorter and without PositiveProjectionsEliminator:
        if (current_mult != 1) {
            subspace.properties.degeneracy /= 2;
        }

        vector_result.emplace_back(std::move(subspace));
    }

    return space::Space(std::move(vector_result));
}

std::unique_ptr<quantum::linear_algebra::AbstractDenseSemiunitaryMatrix>
S2Transformer::constructTransformationMatrix(
    const std::vector<spin_algebra::SSquaredState>& s_squared_states,
    const Subspace& subspace) const {
    const auto number_of_sz_states = subspace.decomposition->size_cols();
    const auto number_of_s2_states = s_squared_states.size();
    auto transformation_matrix =
        factories_.createDenseSemiunitaryMatrix(number_of_s2_states, number_of_sz_states);
    const auto& ssquared_converter = *ssquared_converter_;

#pragma omp parallel for shared( \
            number_of_s2_states, \
            number_of_sz_states, \
            s_squared_states, \
            ssquared_converter, \
            subspace, \
            transformation_matrix) default(none) schedule(static)
    for (size_t j = 0; j < number_of_sz_states; ++j) {
        auto iterator = subspace.decomposition->GetNewIterator(j);
        std::vector<std::vector<double>> all_projections;
        std::vector<double> all_coeffs;
        all_projections.reserve(iterator->size());
        all_coeffs.reserve(iterator->size());

        while (iterator->hasNext()) {
            auto item = iterator->getNext();
            auto lex_index = item.index;

            all_projections.emplace_back(std::move(construct_projections(lex_index)));
            all_coeffs.emplace_back(item.value);
        }

        for (size_t i = 0; i < number_of_s2_states; ++i) {
            const auto& s_squared_state = s_squared_states[i];
            double value = 0;
            double acc = 0;

            // iterate over all lex-states of these (possibly) symmetrized state
            for (size_t k = 0; k < all_coeffs.size(); ++k) {
                const auto& projections = all_projections.at(k);
                const double coeff = all_coeffs.at(k);
                value += ssquared_converter.total_CG_coefficient(s_squared_state, projections) / coeff;
                acc++;
            }
            // something like "averaging" of value:
            value /= acc;
            // I do not understand why, but it works.
            transformation_matrix->add_to_position(value, i, j);
        }
    }

    // Suppose, we have intermediate states : |...01> and |...10>.
    // If some S2-group mixes these states, we actually need to use
    // (|...01> + |...10>) / sqrt(2) and (|...01> - |...10>) / sqrt(2) instead.
    // To avoid explicit using of these states,
    // we will calculate coefficient only for one of states.
    // It breaks norm of semiunitary matrix, so we need to fix it:
    transformation_matrix->normalize();
    // I guess, it works for S2, but I'm not sure, if it will work for non-Abelian groups.

    common::Logger::verbose(
        "Matrix ({0:>6}, {1:>6}) was constructed.",
        number_of_s2_states,
        number_of_sz_states);
    return transformation_matrix;
}

std::vector<double> S2Transformer::construct_projections(uint32_t lex_index) const {
    const auto number_of_mults = converter_->get_mults().size();

    std::vector<double> projections;
    projections.resize(2 * number_of_mults - 1);
    for (size_t a = 0; a < number_of_mults; ++a) {
        projections[a] = converter_->convert_lex_index_to_one_sz_projection(lex_index, a)
            - converter_->get_spins()[a];
    }

    for (const auto& instruction : *ssquared_converter_->getOrderOfSummation()) {
        projections[instruction.position_of_sum] = 0;
        for (const auto& pos : instruction.positions_of_summands) {
            projections[instruction.position_of_sum] += projections[pos];
        }
    }

    return projections;
}

}  // namespace space::optimization