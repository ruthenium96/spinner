#include "MeanAndDeviation.h"

#include <cmath>
#include <iostream>

std::pair<std::chrono::duration<double>, std::chrono::duration<double>>
meanAndDeviation(const std::vector<std::chrono::duration<double>>& cycle_times) {
    double mean_s = 0.0;
    for (const auto& time : cycle_times) {
        mean_s += time.count();
    }
    mean_s /= double(cycle_times.size());
    auto mean = std::chrono::duration<double>(mean_s);

    double standardDeviation_s = 0.0;
    for (const auto& time : cycle_times) {
        double diff_s = std::abs(time.count() - mean.count());
        standardDeviation_s += diff_s * diff_s;
    }
    standardDeviation_s /= double(cycle_times.size());
    standardDeviation_s = sqrt(standardDeviation_s);

    auto standardDeviation = std::chrono::duration<double>(standardDeviation_s);
    return {mean, standardDeviation};
}

void PerformanceTest(const std::function<void(void)>& f, uint32_t cycles) {
    std::vector<std::chrono::duration<double>> cycle_times(cycles);
    for (int i = 0; i < cycles; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        f();
        auto finish = std::chrono::high_resolution_clock::now();
        cycle_times[i] = (finish - start);
    }
    auto [mean, standardDeviation] = meanAndDeviation(cycle_times);
    std::cout << mean.count() << " ± " << standardDeviation.count() << " s\n";
}