#include "Operator.h"

#include "ConstantTerm.h"
#include "ScalarProductTerm.h"

namespace model::operators {
Operator Operator::s_squared(const std::vector<double>& spins) {
    Operator s_squared_operator_;
    double sum_of_s_squared = 0;
    for (double spin : spins) {
        sum_of_s_squared += spin * (spin + 1);
    }
    s_squared_operator_.zero_center_terms.emplace_back(
        std::make_unique<const ConstantTerm>(sum_of_s_squared));
    s_squared_operator_.two_center_terms.emplace_back(
        std::make_unique<const ScalarProductTerm>(spins.size()));
    return s_squared_operator_;
}

Operator::Operator(const Operator& rhs) {
    for (const std::unique_ptr<const ZeroCenterTerm>& el : rhs.zero_center_terms) {
        zero_center_terms.emplace_back(el->clone());
    }
    for (const std::unique_ptr<const OneCenterTerm>& el : rhs.one_center_terms) {
        one_center_terms.emplace_back(el->clone());
    }
    for (const std::unique_ptr<const TwoCenterTerm>& el : rhs.two_center_terms) {
        two_center_terms.emplace_back(el->clone());
    }
}
}  // namespace model::operators