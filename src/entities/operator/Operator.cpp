#include "Operator.h"

#include "src/components/operator/ConstantOperator.h"
#include "src/components/operator/ScalarProduct.h"

Operator Operator::s_squared(const std::vector<double>& spins) {
    Operator s_squared_operator_;
    double sum_of_s_squared = 0;
    for (double spin : spins) {
        sum_of_s_squared += spin * (spin + 1);
    }
    s_squared_operator_.zero_center_terms.emplace_back(
        std::make_unique<const ConstantOperator>(sum_of_s_squared));
    s_squared_operator_.two_center_terms.emplace_back(
        std::make_unique<const ScalarProduct>(spins.size()));
    return s_squared_operator_;
}