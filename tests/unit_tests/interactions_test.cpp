#include <cmath>

#include "gtest/gtest.h"
#include "src/common/lexicographic/IndexConverter.h"
#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/matrix/Matrix.h"
#include "src/model/operators/terms/ConstantTerm.h"
#include "src/model/operators/terms/ScalarProductTerm.h"

TEST(constant_operator, 2222_333_2345_44444) {
    std::vector<std::vector<int>> vector_of_mults =
        {{2, 2, 2, 2}, {3, 3, 3}, {2, 3, 4, 5}, {4, 4, 4, 4, 4}};
    for (const auto& mults : vector_of_mults) {
        // Construct Converter
        lexicographic::IndexConverter converter = lexicographic::IndexConverter(mults);

        auto factories_ = quantum::linear_algebra::FactoriesList();

        // Construct Space
        space::Space space_(converter.get_total_space_size(), factories_);

        for (int constant = 0; constant < 100; constant += 11) {
            auto constant_ptr = std::make_shared<double>(constant);
            // Construct Operator
            model::operators::Operator operator_;
            operator_.getZeroCenterTerms().emplace_back(
                std::make_unique<model::operators::ConstantTerm>(constant_ptr));

            // Build Matrix
            Matrix matrix = Matrix(space_, operator_, converter, factories_);

            // Check results:
            for (const auto& matrix_block : matrix.blocks) {
                for (size_t i = 0; i < matrix_block.raw_data->size(); ++i) {
                    for (size_t j = 0; j < matrix_block.raw_data->size(); ++j) {
                        if (i == j) {
                            EXPECT_DOUBLE_EQ(matrix_block.raw_data->at(i, i), constant);
                        } else {
                            EXPECT_DOUBLE_EQ(matrix_block.raw_data->at(i, j), 0);
                        }
                    }
                }
            }
        }
    }
}

TEST(scalar_product, one_center_1_2_3_4_5_6) {
    // Construct matrix of interaction parameters
    auto ptr_to_js = std::make_shared<const model::TwoDNumericalParameters<double>>(1, NAN);

    std::vector<std::vector<int>> vector_of_mults = {{1}, {2}, {3}, {4}, {5}, {6}};
    for (const auto& mults : vector_of_mults) {
        // Construct Converter
        lexicographic::IndexConverter converter = lexicographic::IndexConverter(mults);

        auto factories_ = quantum::linear_algebra::FactoriesList();

        // Construct Space
        space::Space space_(converter.get_total_space_size(), factories_);

        // Construct Operator
        model::operators::Operator operator_;
        operator_.getTwoCenterTerms().emplace_back(
            std::make_unique<model::operators::ScalarProductTerm>(converter, ptr_to_js));

        // Build Matrix
        Matrix matrix = Matrix(space_, operator_, converter, factories_);

        // Check results:
        for (const auto& matrix_block : matrix.blocks) {
            for (size_t i = 0; i < matrix_block.raw_data->size(); ++i) {
                for (size_t j = 0; j < matrix_block.raw_data->size(); ++j) {
                    EXPECT_DOUBLE_EQ(matrix_block.raw_data->at(i, j), 0);
                }
            }
        }
    }
}

TEST(scalar_product, one_interaction_22_222_2222_33_333_3333_44_444_4444_23456) {
    std::vector<std::vector<int>> vector_of_mults = {
        {2, 2},
        {2, 2, 2},
        {2, 2, 2, 2},
        {3, 3},
        {3, 3, 3},
        {3, 3, 3, 3},
        {4, 4},
        {4, 4, 4},
        {4, 4, 4, 4},
        {2, 3, 4, 5, 6}};
    for (const auto& mults : vector_of_mults) {
        // Construct matrix of interaction parameters
        double J = 10;
        auto ptr_to_js =
            std::make_shared<model::TwoDNumericalParameters<double>>(mults.size(), NAN);
        ptr_to_js->at(0, 1) = J;

        // Construct Converter
        lexicographic::IndexConverter converter = lexicographic::IndexConverter(mults);

        auto factories_ = quantum::linear_algebra::FactoriesList();

        // Construct Space
        space::Space space_(converter.get_total_space_size(), factories_);

        // Construct Operator
        model::operators::Operator operator_;
        operator_.getTwoCenterTerms().emplace_back(
            std::make_unique<model::operators::ScalarProductTerm>(converter, ptr_to_js));

        // Create Factory:
        auto factory_ = quantum::linear_algebra::AbstractSymmetricMatrixFactory::defaultFactory();

        // Build Matrix
        Matrix matrix = Matrix(space_, operator_, converter, factories_);

        // Check results:
        for (const auto& matrix_block : matrix.blocks) {
            for (size_t i = 0; i < matrix_block.raw_data->size(); ++i) {
                for (size_t j = 0; j < matrix_block.raw_data->size(); ++j) {
                    double first_spin = converter.get_spins()[0];
                    double second_spin = converter.get_spins()[1];
                    double j_first_center_projection =
                        converter.convert_lex_index_to_one_sz_projection(j, 0) - first_spin;
                    double j_second_center_projection =
                        converter.convert_lex_index_to_one_sz_projection(j, 1) - second_spin;
                    // NB: ScalarProduct now use -2*J coefficients.
                    if (i == j) {
                        // S0z * S1z part
                        EXPECT_DOUBLE_EQ(
                            matrix_block.raw_data->at(i, i),
                            -2 * J * j_first_center_projection * j_second_center_projection);
                    } else if (
                        i
                        == converter
                               .ladder_projection(converter.ladder_projection(j, 0, +1), 1, -1)) {
                        // S0+ * S1- part
                        // sqrt(S * (S + 1) - M * (M + 1))
                        double first_center_ladder_factor = sqrt(
                            first_spin * (first_spin + 1)
                            - j_first_center_projection * (j_first_center_projection + 1));
                        // sqrt(S * (S + 1) - M * (M - 1))
                        double second_center_ladder_factor = sqrt(
                            second_spin * (second_spin + 1)
                            - j_second_center_projection * (j_second_center_projection - 1));
                        EXPECT_DOUBLE_EQ(
                            matrix_block.raw_data->at(i, j),
                            -J * first_center_ladder_factor * second_center_ladder_factor);
                    } else if (
                        i
                        == converter
                               .ladder_projection(converter.ladder_projection(j, 0, -1), 1, +1)) {
                        // S0- * S1+ part
                        double first_center_ladder_factor = sqrt(
                            first_spin * (first_spin + 1)
                            - j_first_center_projection * (j_first_center_projection - 1));
                        double second_center_ladder_factor = sqrt(
                            second_spin * (second_spin + 1)
                            - j_second_center_projection * (j_second_center_projection + 1));
                        EXPECT_DOUBLE_EQ(
                            matrix_block.raw_data->at(i, j),
                            -J * first_center_ladder_factor * second_center_ladder_factor);
                    } else {
                        EXPECT_DOUBLE_EQ(matrix_block.raw_data->at(i, j), 0);
                    }
                }
            }
        }
    }
}