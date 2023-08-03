#include <cmath>
#include <deque>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"
#include "src/spin_algebra/MultiplicityDirectSum.h"

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
            runner.getOperator(common::QuantityEnum::S_total_squared),
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

// TODO: add tests for OrderOfSummation