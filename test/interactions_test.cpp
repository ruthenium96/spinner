#include "gtest/gtest.h"

#include "common/LexicographicIndexConverter.h"
#include "components/matrix/MatrixBuilder.h"
#include "components/operator/ConstantOperator.h"
#include "components/operator/ScalarProduct.h"
#include "entities/operator/Operator.h"

TEST(constant_operator, 2222_333_2345_44444) {

    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2},
                                                     {3, 3, 3},
                                                     {2, 3, 4, 5},
                                                     {4, 4, 4, 4, 4}};
    for (const auto& mults : vector_of_mults) {
        // Construct Converter
        spaces::LexicographicIndexConverter converter = spaces::LexicographicIndexConverter(mults);

        // Construct Space
        Space space_(converter.total_space_size);

        for (int constant = 0; constant < 100; constant += 11) {
            // Construct Operator
            Operator operator_;
            operator_.zero_center_terms.emplace_back(new ConstantOperator(constant));

            // Build Matrix
            MatrixBuilder matrixBuilder(converter);
            Matrix matrix = matrixBuilder.apply(space_, operator_);

            // Check results:
            for (const auto& matrix_block : matrix.blocks) {
                for (size_t i = 0; i < matrix_block.raw_data.size(); ++i) {
                    for (size_t j = 0; j < matrix_block.raw_data.size(); ++j) {
                        if (i == j) {
                            EXPECT_DOUBLE_EQ(matrix_block.raw_data(i, i), constant);
                        } else {
                            EXPECT_DOUBLE_EQ(matrix_block.raw_data(i, j), 0);
                        }
                    }
                }
            }
        }
    }
}

TEST(scalar_product, one_center_1_2_3_4_5_6) {
    // Construct matrix of interaction parameters
    arma::dmat js(1, 1);
    for (size_t i = 0; i < 1; ++i) {
        for (size_t j = 0; j < 1; ++j) {
            js(i, j) = NAN;
        }
    }

    std::vector<std::vector<int>> vector_of_mults = {{1},
                                                     {2},
                                                     {3},
                                                     {4},
                                                     {5},
                                                     {6}};
    for (const auto& mults : vector_of_mults) {
        // Construct Converter
        spaces::LexicographicIndexConverter converter = spaces::LexicographicIndexConverter(mults);

        // Construct Space
        Space space_(converter.total_space_size);

        // Construct Operator
        Operator operator_;
        operator_.two_center_terms.emplace_back(new ScalarProduct(js));

        // Build Matrix
        MatrixBuilder matrixBuilder(converter);
        Matrix matrix = matrixBuilder.apply(space_, operator_);

        // Check results:
        for (const auto& matrix_block : matrix.blocks) {
            for (size_t i = 0; i < matrix_block.raw_data.size(); ++i) {
                for (size_t j = 0; j < matrix_block.raw_data.size(); ++j) {
                    EXPECT_DOUBLE_EQ(matrix_block.raw_data(i, j), 0);
                }
            }
        }
    }
}

TEST(scalar_product, one_interaction_22_222_2222_33_333_3333_44_444_4444_23456) {

    std::vector<std::vector<int>> vector_of_mults = {{2, 2},
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
        arma::dmat js(mults.size(), mults.size());
        for (size_t i = 0; i < mults.size(); ++i) {
            for (size_t j = 0; j < mults.size(); ++j) {
                js(i, j) = NAN;
            }
        }
        double J = 10;
        js(0, 1) = J; js(1, 0) = J;

        // Construct Converter
        spaces::LexicographicIndexConverter converter = spaces::LexicographicIndexConverter(mults);

        // Construct Space
        Space space_(converter.total_space_size);

        // Construct Operator
        Operator operator_;
        operator_.two_center_terms.emplace_back(new ScalarProduct(js));

        // Build Matrix
        MatrixBuilder matrixBuilder(converter);
        Matrix matrix = matrixBuilder.apply(space_, operator_);

        // Check results:
        for (const auto& matrix_block : matrix.blocks) {
            for (size_t i = 0; i < matrix_block.raw_data.size(); ++i) {
                for (size_t j = 0; j < matrix_block.raw_data.size(); ++j) {
                    double first_spin = converter.get_spins()[0];
                    double second_spin = converter.get_spins()[1];
                    double j_first_center_projection = converter.convert_lex_index_to_one_sz_projection(j, 0) - first_spin;
                    double j_second_center_projection = converter.convert_lex_index_to_one_sz_projection(j, 1) - second_spin;
                    // NB: ScalarProduct now use -2*J coefficients.
                    if (i == j) {
                        // S0z * S1z part
                        EXPECT_DOUBLE_EQ(matrix_block.raw_data(i, i), -2 * J * j_first_center_projection * j_second_center_projection);
                    } else if (i == converter.ladder_projection(converter.ladder_projection(j, 0, +1), 1, -1)) {
                        // S0+ * S1- part
                        // sqrt(S * (S + 1) - M * (M + 1))
                        double first_center_ladder_factor = sqrt(first_spin * (first_spin + 1) -
                                j_first_center_projection * (j_first_center_projection + 1));
                        // sqrt(S * (S + 1) - M * (M - 1))
                        double second_center_ladder_factor = sqrt(second_spin * (second_spin + 1) -
                                j_second_center_projection * (j_second_center_projection - 1));
                        EXPECT_DOUBLE_EQ(matrix_block.raw_data(i, j), -J * first_center_ladder_factor * second_center_ladder_factor);
                    } else if (i == converter.ladder_projection(converter.ladder_projection(j, 0, -1), 1, +1)) {
                        // S0- * S1+ part
                        double first_center_ladder_factor = sqrt(first_spin * (first_spin + 1) -
                                j_first_center_projection * (j_first_center_projection - 1));
                        double second_center_ladder_factor = sqrt(second_spin * (second_spin + 1) -
                                j_second_center_projection * (j_second_center_projection + 1));
                        EXPECT_DOUBLE_EQ(matrix_block.raw_data(i, j), -J * first_center_ladder_factor * second_center_ladder_factor);
                    } else {
                        EXPECT_DOUBLE_EQ(matrix_block.raw_data(i, j), 0);
                    }
                }
            }
        }
    }
}