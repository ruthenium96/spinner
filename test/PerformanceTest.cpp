#include "gtest/gtest.h"
#include "entities/space/Space.h"
#include "components/space/Symmetrizer.h"
#include "components/space/TzSorter.h"
#include "common/Logger.h"
#include <chrono>
#include "common/runner/Runner.h"

void PerformanceTest(std::function<void(void)> f, int cycles = 1) {
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < cycles; ++i) {
        f();
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_once = (finish - start) / cycles;
    std::cout << "Elapsed time: " << elapsed_once.count() << " s\n";
}

TEST(performanceTest, simple2ComponentSchema) {

    std::vector<int> mults = {4, 4, 4, 4, 4, 4, 4, 4, 4,};

    lexicographic::IndexConverter converter(mults);

    PerformanceTest([&mults]() {
        runner::Runner runner(mults);

        runner.TzSort();

        runner.Symmetrize(Group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6,}, {0, 2, 1, 3, 5, 4, 6, 8, 7, }});
        runner.Symmetrize(Group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2, }, {0, 1, 2, 6, 7, 8, 3, 4, 5, }});
        }, 1);
}