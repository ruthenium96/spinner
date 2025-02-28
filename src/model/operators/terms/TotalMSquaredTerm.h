#ifndef TOTALMSQUAREDTERM_H
#define TOTALMSQUAREDTERM_H

#include "Term.h"
#include "src/common/index_converter/AbstractIndexConverter.h"

namespace model::operators {
class TotalMSquaredTerm : public model::operators::ZeroCenterTerm{
public:
    TotalMSquaredTerm(
        std::shared_ptr<const index_converter::AbstractIndexConverter> converter);
    std::unique_ptr<Term> clone() const override;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        const std::set<unsigned int>& indexes_of_vectors) const override;

private:
    std::shared_ptr<const index_converter::AbstractIndexConverter> converter_;

};
} // namespace model::operators

#endif