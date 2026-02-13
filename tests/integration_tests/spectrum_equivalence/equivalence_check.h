#ifndef SPINNER_EQUIVALENCECHECK_H
#define SPINNER_EQUIVALENCECHECK_H

#include <string>
#include <vector>
#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"
#include "src/common/physical_optimization/OptimizationList.h"

void expect_final_vectors_equivalence(spinner::runner::Runner& simple, spinner::runner::Runner& second, bool assert=false);

using TestParam = std::tuple<spinner::common::physical_optimization::OptimizationList, std::vector<spinner::spin_algebra::Multiplicity>>;
class SpectrumFinalEquivalenceTest : public ::testing::TestWithParam<TestParam> {};

namespace spinner::common::physical_optimization {
std::ostream& operator<<(std::ostream& os, const OptimizationList& optimization_list);
} // namespace spinner::common::physical_optimization
std::string spectrum_final_equivalence_test_name_generator(const ::testing::TestParamInfo<TestParam>& info);

#endif  //SPINNER_EQUIVALENCECHECK_H