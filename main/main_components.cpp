#include "common/Space.h"
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

int main() {

    std::vector<int> mults = {2, 2, 2, 2};

    PerformanceTest([&mults]() {
        Space T(mults);
        Tz_Sorter tz_sorter(mults);
        T = tz_sorter(T);
        Symmetrizer c2_symmetrizer(mults, 2);
        T = c2_symmetrizer(T);
        T.print();
    });
}