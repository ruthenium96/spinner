#include "Operator.h"

#include "src/model/operators/terms/ConstantTerm.h"
#include "src/model/operators/terms/ScalarProductTerm.h"
#include "src/model/operators/terms/ZProjectionsProductTerm.h"

namespace model::operators {
Operator Operator::s_squared(lexicographic::IndexConverter converter) {
    Operator s_squared_operator_;
    double sum_of_s_squared = 0;
    for (double spin : converter.get_spins()) {
        sum_of_s_squared += spin * (spin + 1);
    }
    s_squared_operator_.zero_center_terms.emplace_back(
        std::make_unique<const ConstantTerm>(sum_of_s_squared));
    s_squared_operator_.two_center_terms.emplace_back(
        std::make_unique<const ScalarProductTerm>(converter));
    return s_squared_operator_;
}

Operator Operator::g_sz_squared(
    lexicographic::IndexConverter converter,
    std::shared_ptr<const TwoDNumericalParameters<double>> parameters) {
    Operator g_sz_squared_operator_;
    g_sz_squared_operator_.two_center_terms.emplace_back(
        std::make_unique<const ZProjectionsProductTerm>(converter, parameters));
    return g_sz_squared_operator_;
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

bool Operator::empty() const {
    return zero_center_terms.empty() && one_center_terms.empty() && two_center_terms.empty();
}

std::vector<std::unique_ptr<const ZeroCenterTerm>>& Operator::getZeroCenterTerms() {
    return zero_center_terms;
}

const std::vector<std::unique_ptr<const ZeroCenterTerm>>& Operator::getZeroCenterTerms() const {
    return zero_center_terms;
}

std::vector<std::unique_ptr<const OneCenterTerm>>& Operator::getOneCenterTerms() {
    return one_center_terms;
}

const std::vector<std::unique_ptr<const OneCenterTerm>>& Operator::getOneCenterTerms() const {
    return one_center_terms;
}

std::vector<std::unique_ptr<const TwoCenterTerm>>& Operator::getTwoCenterTerms() {
    return two_center_terms;
}

const std::vector<std::unique_ptr<const TwoCenterTerm>>& Operator::getTwoCenterTerms() const {
    return two_center_terms;
}

}  // namespace model::operators