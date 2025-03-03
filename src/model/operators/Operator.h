#ifndef SPINNER_OPERATOR_H
#define SPINNER_OPERATOR_H

#include <memory>
#include <vector>

#include "src/model/operators/terms/Term.h"

namespace model::operators {
class Operator {
  public:
    Operator() = default;
    Operator(const Operator&) = delete;
    Operator& operator=(const Operator&) = delete;
    Operator(Operator&&) noexcept = default;
    Operator& operator=(Operator&&) noexcept = default;
    ~Operator() = default;

    bool empty() const;
    const std::vector<std::unique_ptr<const Term>>& getTerms() const;

    void emplace_back(std::unique_ptr<const Term>&& term);

  private:
    std::vector<std::unique_ptr<const Term>> terms_;
};
}  // namespace model::operators

#endif  //SPINNER_OPERATOR_H
