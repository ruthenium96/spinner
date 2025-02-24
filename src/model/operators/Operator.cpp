#include "Operator.h"

namespace model::operators {

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