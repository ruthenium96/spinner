#ifndef SPINNER_CONSTANTTERM_H
#define SPINNER_CONSTANTTERM_H

#include "Term.h"

namespace model::operators {
class ConstantTerm: public ZeroCenterTerm {
  public:
    explicit ConstantTerm(std::shared_ptr<const double> constant);

    std::unique_ptr<Term> clone() const override;

    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        const std::set<unsigned int>& indexes_of_vectors) const override;

  private:
    std::shared_ptr<const double> constant_;
};
}  // namespace model::operators

#endif  //SPINNER_CONSTANTTERM_H
