#ifndef JULY_CONSTANTOPERATOR_H
#define JULY_CONSTANTOPERATOR_H

#include "components/operator/Interaction.h"

class ConstantOperator : public ZeroCenterTerm {
public:
    explicit ConstantOperator(double constant);

    void construct(LexicograficalMatrix& matrix_in_lexicografical_basis, const spaces::LexicographicIndexConverter& converter,
                   uint32_t index_of_vector) const;

private:
    double constant_;
};


#endif //JULY_CONSTANTOPERATOR_H
