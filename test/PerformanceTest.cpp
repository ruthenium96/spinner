#include "gtest/gtest.h"
#include "entities/Space.h"
#include "components/symmetrizer.h"
#include "components/tz_sorter.h"
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

    std::vector<int> mults = {3, 3, 3};

    spaces::LexicographicIndexConverter converter(mults);

    PerformanceTest([&converter]() {
        Space space(converter.total_space_size);
        Tz_Sorter tz_sorter(converter);
        space = tz_sorter.apply(space);
        Group group_first(group::S3, {{1, 2, 0}, {0, 2, 1}});
        Symmetrizer symmetrizer_first(converter, group_first);
        space = symmetrizer_first.apply(space);
        std::cout << space;
    });
}