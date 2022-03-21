#ifndef JULY_OPERATOR_H
#define JULY_OPERATOR_H

#include <memory>
#include <vector>

#include "src/components/operator/Interaction.h"

struct Operator {
    std::vector<std::unique_ptr<const ZeroCenterTerm>> zero_center_terms;
    std::vector<std::unique_ptr<const OneCenterTerm>> one_center_terms;
    std::vector<std::unique_ptr<const TwoCenterTerm>> two_center_terms;

    static Operator s_squared(const std::vector<double>& spins);

    bool empty() const {
        return zero_center_terms.empty() && one_center_terms.empty() && two_center_terms.empty();
    }
};

#endif  //JULY_OPERATOR_H
