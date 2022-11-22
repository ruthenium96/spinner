#include <chrono>

#include "gtest/gtest.h"
#include "src/common/Logger.h"
#include "src/common/runner/Runner.h"
#include "src/space/Space.h"
#include "src/space/optimization/Symmetrizer.h"
#include "src/space/optimization/TzSorter.h"
#include "tests/tools/MeanAndDeviation.h"

TEST(performanceTest, simple2ComponentSchema) {
    std::vector<int> mults = {
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
    };

    lexicographic::IndexConverter converter(mults);

    PerformanceTest(
        [&mults]() {
            model::ModelInput model(mults);
            common::physical_optimization::OptimizationList optimizationList;
            optimizationList
                .Symmetrize(
                    group::Group::S3,
                    {{
                         1,
                         2,
                         0,
                         4,
                         5,
                         3,
                         7,
                         8,
                         6,
                     },
                     {
                         0,
                         2,
                         1,
                         3,
                         5,
                         4,
                         6,
                         8,
                         7,
                     }})
                .Symmetrize(
                    group::Group::S3,
                    {{
                         3,
                         4,
                         5,
                         6,
                         7,
                         8,
                         0,
                         1,
                         2,
                     },
                     {
                         0,
                         1,
                         2,
                         6,
                         7,
                         8,
                         3,
                         4,
                         5,
                     }});
            runner::Runner runner(model, optimizationList);
        },
        1);
}