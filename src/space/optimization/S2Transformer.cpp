#include "S2Transformer.h"

#include <omp.h>

#include <utility>

#include "wignerSymbols.h"

namespace space::optimization {

S2Transformer::S2Transformer(
    lexicographic::IndexConverter converter,
    quantum::linear_algebra::FactoriesList factories,
    std::shared_ptr<const spin_algebra::OrderOfSummation> order_of_summation,
    const spin_algebra::RepresentationsMultiplier& representationsMultiplier) :
    converter_(std::move(converter)),
    factories_(std::move(factories)),
    order_of_summation_(std::move(order_of_summation)) {
    sorted_s_squared_states_ = spin_algebra::SSquaredState::addAllMultiplicitiesAndSort(
        converter_.get_mults(),
        order_of_summation_,
        representationsMultiplier);
    auto max_number_of_threads = omp_get_max_threads();
    hashed_CGs_for_all_threads =
        std::vector<std::unordered_map<uint64_t, double>>(max_number_of_threads);
}

space::Space S2Transformer::apply(Space&& space) const {
    std::vector<Subspace> vector_result;
    vector_result.reserve(space.getBlocks().size());

    std::cout << "S2-transformation started." << std::endl;

    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        Subspace& subspace = space.getBlocks()[i];

        // It is total_proj + max_total_spin = ntz_value
        // We need 2 * |total_proj| + 1 = 2 * |ntz_value - max_total_spin| + 1 =
        // |2 * ntz_value - 2 * max_total_spin - 1 + 1| + 1 =
        // |2 * ntz_value - max_total_proj + 1| + 1
        auto ntz_value = subspace.properties.n_proj.value();
        int a = 2 * (int)ntz_value + 1 - (int)converter_.get_max_ntz_proj();
        spin_algebra::Multiplicity current_mult = std::abs(a) + 1;

        spin_algebra::SSquaredState::Properties subspace_properties;
        subspace_properties.multiplicity = current_mult;
        subspace_properties.representations = subspace.properties.representation;

        if (sorted_s_squared_states_.count(subspace_properties) == 0) {
            continue;
        }
        const auto& s_squared_states = sorted_s_squared_states_.at(subspace_properties);

        subspace.dense_semiunitary_matrix =
            constructTransformationMatrix(s_squared_states, subspace);

        subspace.properties.total_mult = current_mult;
        subspace.properties.degeneracy *= current_mult;
        // TODO: for cases with TzSorter and without PositiveProjectionsEliminator:
        if (current_mult != 1) {
            subspace.properties.degeneracy /= 2;
        }

        vector_result.emplace_back(std::move(subspace));
    }

    std::cout << "S2-transformation finished." << std::endl;

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

#pragma omp parallel for shared( \
        number_of_s2_states, \
            number_of_sz_states, \
            s_squared_states, \
            subspace, \
            transformation_matrix, \
            hashed_CGs_for_all_threads) default(none)
    for (size_t j = 0; j < number_of_sz_states; ++j) {
        auto iterator = subspace.decomposition->GetNewIterator(j);
        // todo: consider reserve or resize
        std::vector<std::vector<double>> all_projections;
        std::vector<double> all_coeffs;

        auto number_of_current_thread = omp_get_thread_num();
        std::unordered_map<uint64_t, double>& hashed_CGs_per_thread =
            hashed_CGs_for_all_threads.at(number_of_current_thread);

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
                value += total_CG_coefficient(s_squared_state, projections, hashed_CGs_per_thread)
                    / coeff;
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

    std::cout << "S2-transformation matrix (" << number_of_s2_states << ", " << number_of_sz_states
              << ") was successfully constructed." << std::endl;
    return transformation_matrix;
}

double S2Transformer::total_CG_coefficient(
    const spin_algebra::SSquaredState& s_squared_state,
    const std::vector<double>& projections,
    std::unordered_map<uint64_t, double>& cached_CGs) const {
    double c = 1;

    for (const auto& instruction : *order_of_summation_) {
        size_t pos_one = instruction.positions_of_summands[0];
        size_t pos_two = instruction.positions_of_summands[1];
        size_t pos_sum = instruction.position_of_sum;

        double spin_one = ((double)s_squared_state.getMultiplicity(pos_one) - 1.0) / 2.0;
        double spin_two = ((double)s_squared_state.getMultiplicity(pos_two) - 1.0) / 2.0;
        double spin_sum = ((double)s_squared_state.getMultiplicity(pos_sum) - 1.0) / 2.0;

        double proj_one = projections[pos_one];
        double proj_two = projections[pos_two];
        c *= hashed_clebsh_gordan(spin_one, spin_two, spin_sum, proj_one, proj_two, cached_CGs);
        if (c == 0.0) {
            return c;
        }
    }

    return c;
}

double S2Transformer::hashed_clebsh_gordan(
    double l1,
    double l2,
    double l3,
    double m1,
    double m2,
    std::unordered_map<uint64_t, double>& cached_CGs) const {
    if (l1 + l2 < l3 || std::abs(l1 - l2) > l3) {
        return 0;
    }
    auto key = to_key(l1, l2, l3, m1, m2);
    if (cached_CGs.contains(key)) {
        return cached_CGs.at(key);
    } else {
        double value = WignerSymbols::clebschGordan(l1, l2, l3, m1, m2, m1 + m2);
        cached_CGs[key] = value;
        return value;
    }
}

std::vector<double> S2Transformer::construct_projections(uint32_t lex_index) const {
    const auto number_of_mults = converter_.get_mults().size();

    std::vector<double> projections;
    projections.resize(2 * number_of_mults - 1);
    for (size_t a = 0; a < number_of_mults; ++a) {
        projections[a] = converter_.convert_lex_index_to_one_sz_projection(lex_index, a)
            - converter_.get_spins()[a];
    }

    for (const auto& instruction : *order_of_summation_) {
        projections[instruction.position_of_sum] = 0;
        for (const auto& pos : instruction.positions_of_summands) {
            projections[instruction.position_of_sum] += projections[pos];
        }
    }

    return projections;
}

uint64_t S2Transformer::to_key(double l1, double l2, double l3, double m1, double m2) noexcept {
    const uint64_t BASE = 4096;
    auto ul1 = (uint64_t)(2 * l1);
    auto ul2 = (uint64_t)(2 * l2);
    auto ul3 = (uint64_t)(2 * l3);
    auto um1 = (uint64_t)(m1 + l1);
    auto um2 = (uint64_t)(m2 + l2);
    uint64_t answer = 0;
    answer |= ul1;
    answer |= ul2 * BASE;
    answer |= ul3 * BASE * BASE;
    answer |= um1 * BASE * BASE * BASE;
    answer |= um2 * BASE * BASE * BASE * BASE;
    return answer;
}
}  // namespace space::optimization