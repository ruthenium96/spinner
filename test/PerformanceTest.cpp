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

    std::vector<int> mults = {2, 2, 2, 2};

    spaces::LexicographicIndexWorker converter(mults);

    PerformanceTest([&converter]() {
        Space space(converter);
        Tz_Sorter tz_sorter(converter);
        space = tz_sorter(space);
        Group group_first(P2, {{1, 0, 3, 2}});
        Symmetrizer symmetrizer_first(converter, group_first);
        space = symmetrizer_first(space);
        Group group_second(P2, {{3, 2, 1, 0}});
        Symmetrizer symmetrizer_second(converter, group_second);
        space = symmetrizer_second(space);
        std::cout << space;
    });
}