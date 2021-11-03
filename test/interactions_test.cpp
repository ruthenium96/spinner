#include "gtest/gtest.h"

#include "common/LexicographicIndexConverter.h"
#include "components/matrix/MatrixBuilder.h"
#include "components/operator/ConstantOperator.h"
#include "entities/operator/Operator.h"

TEST(constant_operator, 2222_333_2345_44444) {

    // Construct Converter
    std::vector<std::vector<int>> vector_of_mults = {{2, 2, 2, 2},
                                                     {3, 3, 3},
                                                     {2, 3, 4, 5},
                                                     {4, 4, 4, 4, 4}};
    for (const auto& mults : vector_of_mults) {
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
                            EXPECT_EQ(matrix_block.raw_data(i, j), constant);
                        } else {
                            EXPECT_EQ(matrix_block.raw_data(i, j), 0);
                        }
                    }
                }
            }
        }
    }
}
