#ifndef SPINNER_OPERATOR_H
#define SPINNER_OPERATOR_H

#include <memory>
#include <vector>

#include "src/model/operators/terms/Term.h"

namespace model::operators {
class Operator {
  public:
    Operator() = default;
    Operator(const Operator&);
    Operator& operator=(const Operator&) = delete;
    Operator(Operator&&) noexcept = default;
    Operator& operator=(Operator&&) noexcept = default;
    ~Operator() = default;

    static Operator s_squared(lexicographic::IndexConverter converter);
    static Operator g_sz_squared(
        lexicographic::IndexConverter converter,
        std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
        std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters);

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
