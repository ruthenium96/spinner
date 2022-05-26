#ifndef SPINNER_OPERATOR_H
#define SPINNER_OPERATOR_H

#include <memory>
#include <vector>

#include "Term.h"

namespace model::operators {
class Operator {
  public:
    Operator() = default;
    Operator(const Operator&);
    // TODO: should we implement copy assignment?
    Operator& operator=(const Operator&) = delete;
    Operator(Operator&&) noexcept = default;
    Operator& operator=(Operator&&) noexcept = default;
    ~Operator() = default;

    static Operator s_squared(const std::vector<double>& spins);

    bool empty() const;
    std::vector<std::unique_ptr<const ZeroCenterTerm>>& getZeroCenterTerms();
    const std::vector<std::unique_ptr<const ZeroCenterTerm>>& getZeroCenterTerms() const;
    std::vector<std::unique_ptr<const OneCenterTerm>>& getOneCenterTerms();
    const std::vector<std::unique_ptr<const OneCenterTerm>>& getOneCenterTerms() const;
    std::vector<std::unique_ptr<const TwoCenterTerm>>& getTwoCenterTerms();
    const std::vector<std::unique_ptr<const TwoCenterTerm>>& getTwoCenterTerms() const;

  private:
    std::vector<std::unique_ptr<const ZeroCenterTerm>> zero_center_terms;
    std::vector<std::unique_ptr<const OneCenterTerm>> one_center_terms;
    std::vector<std::unique_ptr<const TwoCenterTerm>> two_center_terms;
};
}  // namespace model::operators

#endif  //SPINNER_OPERATOR_H
