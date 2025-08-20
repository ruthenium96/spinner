#ifndef SPINNER_EQUIVALENCECHECK_H
#define SPINNER_EQUIVALENCECHECK_H

#include <string>
#include <vector>
#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"
#include "src/common/physical_optimization/OptimizationList.h"

void expect_final_vectors_equivalence(runner::Runner& simple, runner::Runner& second);

using TestParam = std::tuple<common::physical_optimization::OptimizationList, std::vector<spin_algebra::Multiplicity>>;
class SpectrumFinalEquivalenceTest : public ::testing::TestWithParam<TestParam> {};


namespace common::physical_optimization {
std::ostream& operator<<(std::ostream& os, const common::physical_optimization::OptimizationList& optimization_list);
}
std::string spectrum_final_equivalence_test_name_generator(const ::testing::TestParamInfo<TestParam>& info);

#endif  //SPINNER_EQUIVALENCECHECK_H