#include "Operator.h"

#include "src/model/operators/terms/ConstantTerm.h"
#include "src/model/operators/terms/lexicographic/ScalarProductTerm.h"
#include "src/model/operators/terms/lexicographic/SzSzOneCenterTerm.h"
#include "src/model/operators/terms/lexicographic/SzSzTwoCenterTerm.h"

namespace model::operators {
Operator Operator::s_squared(std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter) {
    Operator s_squared_operator_;
    auto sum_of_s_squared = std::make_shared<double>(0);
    for (double spin : converter->get_spins()) {
        *sum_of_s_squared += spin * (spin + 1);
    }
    s_squared_operator_.terms_.emplace_back(
        std::make_unique<const ConstantTerm>(sum_of_s_squared));
    s_squared_operator_.terms_.emplace_back(
        std::make_unique<const lexicographic::ScalarProductTerm>(converter));
    return s_squared_operator_;
}

Operator Operator::g_sz_squared(
    std::shared_ptr<const index_converter::lexicographic::IndexConverter> converter,
    std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
    std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) {
    Operator g_sz_squared_operator_;
    g_sz_squared_operator_.terms_.emplace_back(
        std::make_unique<const lexicographic::SzSzOneCenterTerm>(converter, diagonal_parameters));
    g_sz_squared_operator_.terms_.emplace_back(
        std::make_unique<const lexicographic::SzSzTwoCenterTerm>(converter, nondiagonal_parameters, 2));
    // this two from summation in Submatrix: \sum_{a=1}^N \sum_{b=a+1}^N
    return g_sz_squared_operator_;
}

bool Operator::empty() const {
    return terms_.empty();
}

const std::vector<std::unique_ptr<const Term>>& Operator::getTerms() const {
    return terms_;
}

void Operator::emplace_back(std::unique_ptr<const Term>&& term) {
    terms_.emplace_back(std::move(term));
}

}  // namespace model::operators