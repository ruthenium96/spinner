#ifndef SPINNER_MEANANDDEVIATION_H
#define SPINNER_MEANANDDEVIATION_H

#include <chrono>
#include <functional>

std::pair<std::chrono::duration<double>, std::chrono::duration<double>>
meanAndDeviation(const std::vector<std::chrono::duration<double>>& cycle_times);

void PerformanceTest(const std::function<void(void)>& f, uint32_t cycles = 1);

#endif  //SPINNER_MEANANDDEVIATION_H
