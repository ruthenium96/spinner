#ifndef JULY_OPERATOR_H
#define JULY_OPERATOR_H

#include <memory>
#include <vector>

#include "src/components/operator/Interaction.h"

struct Operator {
    Operator() = default;
    Operator(const Operator&);
    // TODO: should we implement copy assignment?
    Operator& operator=(const Operator&) = delete;
    Operator(Operator&&) noexcept = default;
    Operator& operator=(Operator&&) noexcept = default;
    ~Operator() = default;

    std::vector<std::unique_ptr<const ZeroCenterTerm>> zero_center_terms;
    std::vector<std::unique_ptr<const OneCenterTerm>> one_center_terms;
    std::vector<std::unique_ptr<const TwoCenterTerm>> two_center_terms;

    static Operator s_squared(const std::vector<double>& spins);

    bool empty() const {
        return zero_center_terms.empty() && one_center_terms.empty() && two_center_terms.empty();
    }
};

#endif  //JULY_OPERATOR_H
