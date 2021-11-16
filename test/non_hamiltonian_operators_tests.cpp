#include <deque>

#include "common/runner/Runner.h"
#include "components/matrix/MatrixBuilder.h"
#include "gtest/gtest.h"

std::deque<int> spin_addition(const std::vector<int>& mults) {
    std::deque<int> total_multiplicities;

    total_multiplicities.push_back(1);

    for (auto mult_of_center : mults) {
        size_t old_size = total_multiplicities.size();
        for (size_t i = 0; i < old_size; ++i) {
            int current_multiplicity = total_multiplicities.front();
            for (int new_multiplicity = std::abs(current_multiplicity - mult_of_center) + 1;
                 new_multiplicity < current_multiplicity + mult_of_center;
                 new_multiplicity += 2) {
                total_multiplicities.push_back(new_multiplicity);
            }
            total_multiplicities.pop_front();
        }
    }
    std::stable_sort(total_multiplicities.begin(), total_multiplicities.end());
    return total_multiplicities;
}

std::vector<int> duplicate_multiplicity_multiplicity_times(const std::deque<int>& multiplicities) {
    std::vector<int> duplicated_multiplicities;
    for (auto mult : multiplicities) {
        for (size_t j = 0; j < mult; ++j) {
            duplicated_multiplicities.push_back(mult);
        }
    }
    std::stable_sort(duplicated_multiplicities.begin(), duplicated_multiplicities.end());
    return duplicated_multiplicities;
}

TEST(spin_addition, 22) {
    std::vector<int> mults = {2, 2};
    std::deque<int> expected_total_multiplicities = {1, 3};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 222) {
    std::vector<int> mults = {2, 2, 2};
    std::deque<int> expected_total_multiplicities = {2, 2, 4};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 2222) {
    std::vector<int> mults = {2, 2, 2, 2};
    std::deque<int> expected_total_multiplicities = {1, 1, 3, 3, 3, 5};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 33) {
    std::vector<int> mults = {3, 3};
    std::deque<int> expected_total_multiplicities = {1, 3, 5};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(spin_addition, 44) {
    std::vector<int> mults = {4, 4};
    std::deque<int> expected_total_multiplicities = {1, 3, 5, 7};
    EXPECT_EQ(spin_addition(mults), expected_total_multiplicities);
}

TEST(
    initialize_s_squared,
    eigenvalues_of_s_squared_matrix_correspond_to_spin_addition_1_22_333_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults =
        {{1}, {2, 2}, {3, 3, 3}, {4, 4, 4, 4}, {2, 3, 4, 5, 6}};

    for (const auto& mults : vector_of_mults) {
        runner::Runner runner(mults);
        lexicographic::IndexConverter converter(mults);

        runner.InitializeSSquared();

        MatrixBuilder matrix_builder(runner.getIndexConverter());
        Matrix s_squared_matrix = matrix_builder.apply(
            runner.getSpace(),
            runner.getOperator(common::QuantityEnum::S_total_squared));

        std::vector<DenseVector> s_squared_values(s_squared_matrix.blocks.size());
        for (size_t i = 0; i < s_squared_values.size(); ++i) {
            s_squared_matrix.blocks[i].raw_data.diagonalize(s_squared_values[i]);
        }
        std::vector<double> s_squared_vector = concatenate(s_squared_values);

        std::vector<int> total_multiplicities(s_squared_vector.size());
        for (size_t i = 0; i < total_multiplicities.size(); ++i) {
            total_multiplicities[i] = (int)round(sqrt(1 + 4 * s_squared_vector[i]));
        }
        std::sort(total_multiplicities.begin(), total_multiplicities.end());

        EXPECT_EQ(
            total_multiplicities,
            duplicate_multiplicity_multiplicity_times(spin_addition(mults)));
    }
}
