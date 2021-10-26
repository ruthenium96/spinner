#include "ConstantOperator.h"

void ConstantOperator::construct(LexicograficalMatrix &matrix_in_lexicografical_basis,
                         const spaces::LexicographicIndexConverter &converter, uint32_t index_of_vector) const {
    matrix_in_lexicografical_basis(index_of_vector, index_of_vector) += constant_;
}

ConstantOperator::ConstantOperator(double constant) : constant_(constant) {
}
