#include <cmath>
#include <deque>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"
#include "src/spin_algebra/GroupAdapter.h"
#include "src/spin_algebra/MultiplicityDirectSum.h"
#include "src/spin_algebra/OrderOfSummation.h"
#include "src/spin_algebra/SSquaredState.h"

spin_algebra::MultiplicityDirectSum
spin_addition(const std::vector<spin_algebra::MultiplicityDirectSum>& mults) {
    return std::accumulate(
        mults.begin(),
        mults.end(),
        spin_algebra::MultiplicityDirectSum(1),  // Identity element of spin addition: S=0.
        std::multiplies<>());
}

std::vector<uint16_t> duplicate_multiplicity_multiplicity_times(
    const spin_algebra::MultiplicityDirectSum& multiplicities) {
    std::vector<uint16_t> duplicated_multiplicities;
    for (auto mult : multiplicities.getMultiplicities()) {
        for (size_t j = 0; j < mult; ++j) {
            duplicated_multiplicities.push_back(mult);
        }
    }
    std::stable_sort(duplicated_multiplicities.begin(), duplicated_multiplicities.end());
    return duplicated_multiplicities;
}

TEST(spin_addition, 22) {
    std::vector<spin_algebra::MultiplicityDirectSum> mults = {{2}, {2}};
    spin_algebra::MultiplicityDirectSum expected_total_multiplicities = {1, 3};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 222) {
    std::vector<spin_algebra::MultiplicityDirectSum> mults = {{2}, {2}, {2}};
    spin_algebra::MultiplicityDirectSum expected_total_multiplicities = {2, 2, 4};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 2222) {
    std::vector<spin_algebra::MultiplicityDirectSum> mults = {{2}, {2}, {2}, {2}};
    spin_algebra::MultiplicityDirectSum expected_total_multiplicities = {1, 1, 3, 3, 3, 5};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 23) {
    std::vector<spin_algebra::MultiplicityDirectSum> mults = {{2}, {3}};
    spin_algebra::MultiplicityDirectSum expected_total_multiplicities = {2, 4};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 33) {
    std::vector<spin_algebra::MultiplicityDirectSum> mults = {{3}, {3}};
    spin_algebra::MultiplicityDirectSum expected_total_multiplicities = {1, 3, 5};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 44) {
    std::vector<spin_algebra::MultiplicityDirectSum> mults = {{4}, {4}};
    spin_algebra::MultiplicityDirectSum expected_total_multiplicities = {1, 3, 5, 7};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

void allMultiplicitiesWereSummedExactlyOnce(
    const spin_algebra::OrderOfSummation& order_of_summation,
    size_t number_of_all_mults) {
    std::vector<bool> multiplicity_was_summed(number_of_all_mults - 1, false);
    for (const auto& instruction : order_of_summation) {
        for (const auto& position : instruction.positions_of_summands) {
            EXPECT_FALSE(multiplicity_was_summed[position])
                << "Multiplicity " << position << " has been summed multiple times.";
            multiplicity_was_summed[position] = true;
        }
    }
    EXPECT_TRUE(std::all_of(
        multiplicity_was_summed.begin(),
        multiplicity_was_summed.end(),
        [](bool b) { return b; }))
        << "Not all multiplicities were summed";
}

void allMultiplicitiesAreReachableLeftToRight(
    const spin_algebra::OrderOfSummation& order_of_summation,
    size_t number_of_initial_mults,
    size_t number_of_all_mults) {
    std::vector<bool> multiplicity_are_reachable(number_of_all_mults, false);
    for (size_t i = 0; i < number_of_initial_mults; ++i) {
        multiplicity_are_reachable[i] = true;
    }

    for (const auto& instruction : order_of_summation) {
        for (const auto& position : instruction.positions_of_summands) {
            EXPECT_TRUE(multiplicity_are_reachable[position])
                << "Multiplicity " << position << " is not reachable.";
        }
        multiplicity_are_reachable[instruction.position_of_sum] = true;
    }
}

std::string representationName(const std::vector<uint8_t>& representations) {
    std::string representation_name;
    for (auto representation : representations) {
        representation_name += representation;
    }
    return representation_name;
}

TEST(order_of_summation, AAAA_S1_S2_S2xS2) {
    group::Group group_one(group::Group::S2, {{1, 0, 3, 2}});
    group::Group group_two(group::Group::S2, {{2, 3, 0, 1}});
    size_t number_of_initial_mults = 4;
    size_t number_of_summation = 3;  // derived it from group information somehow
    size_t number_of_all_mults = number_of_initial_mults + number_of_summation;

    std::vector<std::vector<group::Group>> set_of_groups =
        {{}, {group_one}, {group_two}, {group_one, group_two}, {group_two, group_one}};

    for (const auto& groups : set_of_groups) {
        auto group_adapter = spin_algebra::GroupAdapter(groups, number_of_initial_mults);
        allMultiplicitiesWereSummedExactlyOnce(
            *group_adapter.getOrderOfSummations(),
            number_of_all_mults);
        allMultiplicitiesAreReachableLeftToRight(
            *group_adapter.getOrderOfSummations(),
            number_of_initial_mults,
            number_of_all_mults);
    }
}

TEST(order_of_summation, AAAAAAAAA_S1_S2_S2xS2) {
    group::Group group_one(group::Group::S2, {{6, 7, 8, 3, 4, 5, 0, 1, 2}});
    group::Group group_two(group::Group::S2, {{2, 1, 0, 5, 4, 3, 8, 7, 6}});
    size_t number_of_initial_mults = 9;
    size_t number_of_summation = 8;  // derived it from group information somehow
    size_t number_of_all_mults = number_of_initial_mults + number_of_summation;

    std::vector<std::vector<group::Group>> set_of_groups =
        {{}, {group_one}, {group_two}, {group_one, group_two}, {group_two, group_one}};

    for (const auto& groups : set_of_groups) {
        auto group_adapter = spin_algebra::GroupAdapter(groups, number_of_initial_mults);
        allMultiplicitiesWereSummedExactlyOnce(
            *group_adapter.getOrderOfSummations(),
            number_of_all_mults);
        allMultiplicitiesAreReachableLeftToRight(
            *group_adapter.getOrderOfSummations(),
            number_of_initial_mults,
            number_of_all_mults);
    }
}

// TODO: add more cases and tests for OrderOfSummation

TEST(addAllMultiplicitiesAndSort, 2222) {
    std::vector<spin_algebra::Multiplicity> mults = {2, 2, 2, 2};
    size_t number_of_initial_mults = mults.size();
    size_t number_of_summation = mults.size() - 1;
    group::Group group_one(group::Group::S2, {{1, 0, 3, 2}});
    group::Group group_two(group::Group::S2, {{2, 3, 0, 1}});

    std::vector<std::pair<
        std::vector<group::Group>,
        std::map<spin_algebra::SSquaredState::Properties, size_t>>>
        input_and_correct_answers;

    input_and_correct_answers.push_back({{}, {{{5, {}}, 1}, {{3, {}}, 3}, {{1, {}}, 2}}});

    input_and_correct_answers.push_back(
        {{group_one}, {{{5, {0}}, 1}, {{3, {0}}, 1}, {{3, {1}}, 2}, {{1, {0}}, 2}}});

    input_and_correct_answers.push_back(
        {{group_two}, {{{5, {0}}, 1}, {{3, {0}}, 1}, {{3, {1}}, 2}, {{1, {0}}, 2}}});

    input_and_correct_answers.push_back(
        {{group_one, group_two},
         {{{5, {0, 0}}, 1},
          {{3, {0, 1}}, 1},
          {{3, {1, 0}}, 1},
          {{3, {1, 1}}, 1},
          {{1, {0, 0}}, 2}}});

    for (const auto& [groups, correct_answer] : input_and_correct_answers) {
        auto group_adapter = spin_algebra::GroupAdapter(groups, number_of_initial_mults);
        auto sorted_sum = spin_algebra::SSquaredState::addAllMultiplicitiesAndSort(
            mults,
            group_adapter.getOrderOfSummations(),
            group_adapter.getRepresentationMultiplier());

        for (const auto& [properties, number] : correct_answer) {
            EXPECT_EQ(number, sorted_sum.at(properties).size())
                << "The number of states with multiplicity = " << properties.multiplicity
                << " and representation " << representationName(properties.representations)
                << " should be equal to " << number << ", but actually is equal to "
                << sorted_sum.at(properties).size();
        }
    }
}

TEST(
    initialize_s_squared,
    eigenvalues_of_s_squared_matrix_correspond_to_spin_addition_1_22_333_4444_23456) {
    std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults =
        {{1}, {2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        model::ModelInput model(mults);

        runner::Runner runner(model);
        lexicographic::IndexConverter converter(mults);

        Matrix s_squared_matrix = Matrix(
            runner.getSpace(),
            *runner.getOperator(common::QuantityEnum::S_total_squared).value(),
            runner.getIndexConverter(),
            runner.getDataStructuresFactories());

        auto s_squared_vector = runner.getDataStructuresFactories().createVector();

        for (auto& block : s_squared_matrix.blocks) {
            s_squared_vector->concatenate_with(block.raw_data->diagonalizeValues());
        }

        std::vector<uint16_t> total_multiplicities(s_squared_vector->size());
        for (size_t i = 0; i < total_multiplicities.size(); ++i) {
            total_multiplicities[i] = (uint16_t)round(sqrt(1 + 4 * s_squared_vector->at(i)));
        }
        std::sort(total_multiplicities.begin(), total_multiplicities.end());

        // TODO: refactor it!
        std::vector<spin_algebra::MultiplicityDirectSum> copied_mults;
        copied_mults.reserve(mults.size());
        for (const auto& i : mults) {
            copied_mults.emplace_back(i);
        }

        EXPECT_EQ(
            total_multiplicities,
            duplicate_multiplicity_multiplicity_times(spin_addition(copied_mults)));
    }
}
