#include <chrono>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"
#include "src/space/Space.h"
#include "src/space/optimization/Symmetrizer.h"
#include "tests/tools/AllSparseSemiunitaryMatrixFactories.h"
#include "tests/tools/MeanAndDeviation.h"

TEST(performanceTest, simple2ComponentSchema) {
    std::vector<spin_algebra::Multiplicity> mults = {4, 4, 4, 4, 4, 4, 4, 4, 4};

    lexicographic::IndexConverter converter(mults);

    auto sparseSemiunitaryfactories = constructAllSparseSemiunitaryMatrixFactories();
    std::vector<quantum::linear_algebra::FactoriesList> allFactoriesLists;

    for (const auto& sparseSemiunitaryfactory : sparseSemiunitaryfactories) {
        quantum::linear_algebra::FactoriesList factoryList = quantum::linear_algebra::FactoriesList(
            quantum::linear_algebra::AbstractSymmetricMatrixFactory::defaultFactory(),
            sparseSemiunitaryfactory);
        allFactoriesLists.push_back(factoryList);
    }

    for (const auto& factoryList : allFactoriesLists)
        PerformanceTest(
            [&mults, &factoryList]() {
                model::ModelInput model(mults);
                common::physical_optimization::OptimizationList optimizationList;
                optimizationList.TzSort()
                    .Symmetrize(
                        group::Group::S3,
                        {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}})
                    .Symmetrize(
                        group::Group::S3,
                        {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
                runner::Runner runner(model, optimizationList, factoryList);
            },
        1);
}