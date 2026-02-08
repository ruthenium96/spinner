#include <stdexcept>
#include "gtest/gtest.h"
#include "src/common/physical_optimization/OptimizationList.h"

using namespace spinner::common::physical_optimization;

TEST(ftlm_tests, throw_no_positive_projections_eliminationssquared_transformation) {
    OptimizationList optimizationList;
    optimizationList.TzSort();
    EXPECT_THROW(optimizationList.FTLMApproximate({400, 1024, 200}),
        std::invalid_argument);
}

TEST(ftlm_tests, throw_krylov_size_bigger_than_threshold) {
    OptimizationList optimizationList;
    optimizationList.TzSort().EliminatePositiveProjections();
    EXPECT_THROW(optimizationList.FTLMApproximate({400, 256, 200}),
        std::invalid_argument);
}

TEST(ftlm_tests, throw_number_of_seed_equals_to_zero) {
    OptimizationList optimizationList;
    optimizationList.TzSort().EliminatePositiveProjections();
    EXPECT_THROW(optimizationList.FTLMApproximate({400, 1024, 0}),
        std::invalid_argument);
}