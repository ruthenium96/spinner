#ifndef JULY_OPERATOR_H
#define JULY_OPERATOR_H

#include <vector>

#include "components/operator/Interaction.h"

class Operator {
public:
    std::vector<std::unique_ptr<const ZeroCenterTerm>> zero_center_terms;
    std::vector<std::unique_ptr<const OneCenterTerm>> one_center_terms;
    std::vector<std::unique_ptr<const TwoCenterTerm>> two_center_terms;
private:
};


#endif //JULY_OPERATOR_H
