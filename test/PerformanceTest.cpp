#include "gtest/gtest.h"
#include "entities/Space.h"
#include "components/Symmetrizer.h"
#include "components/TzSorter.h"
#include "common/Logger.h"
#include <chrono>

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

    spaces::LexicographicIndexConverter converter(mults);

    PerformanceTest([&converter]() {

        Space space(converter.total_space_size);

        TzSorter tz_sorter(converter);

        space = tz_sorter.apply(space);

        Group group_S3_first(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6,}, {0, 2, 1, 3, 5, 4, 6, 8, 7, }});
        Symmetrizer symmetrizer_S3_first(converter, group_S3_first);
        space = symmetrizer_S3_first.apply(space);

        Group group_S3_second(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2, }, {0, 1, 2, 6, 7, 8, 3, 4, 5, }});
        Symmetrizer symmetrizer_S3_second(converter, group_S3_second);
        space = symmetrizer_S3_second.apply(space);
        }, 1);
}