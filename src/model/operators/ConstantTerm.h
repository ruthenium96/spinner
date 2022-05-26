#ifndef JULY_CONSTANTTERM_H
#define JULY_CONSTANTTERM_H

#include "Term.h"

namespace model::operators {
class ConstantTerm: public ZeroCenterTerm {
  public:
    explicit ConstantTerm(double constant);

    std::unique_ptr<ZeroCenterTerm> clone() const override;

    void construct(
        lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector) const override;

  private:
    double constant_;
};
}  // namespace model::operators

#endif  //JULY_CONSTANTTERM_H
