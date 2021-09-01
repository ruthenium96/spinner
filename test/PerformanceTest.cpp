#include "gtest/gtest.h"
#include "entities/Space.h"
#include "components/c2_symmetrizer.h"
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

    spaces::LexicographicIndexConverter converter(mults);

    PerformanceTest([&converter]() {
        Space space(converter);
        Tz_Sorter tz_sorter(converter);
        space = tz_sorter(space);
        Symmetrizer c2_symmetrizer(converter, 2);
        space = c2_symmetrizer(space);
        std::cout << space;
    });
}