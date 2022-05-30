#ifndef SPINNER_TERM_H
#define SPINNER_TERM_H

#include <cstdint>
#include <memory>

#include "src/common/lexicographic/IndexConverter.h"
#include "src/common/lexicographic/SparseMatrix.h"
#include "src/entities/data_structures/DenseMatrix.h"

namespace model::operators {
class ZeroCenterTerm {
  public:
    // This method is required for deep copy of std::vector<std::unique_ptr<ZeroCenterTerm>>
    virtual std::unique_ptr<ZeroCenterTerm> clone() const = 0;
    virtual void construct(
        lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector) const = 0;
    virtual ~ZeroCenterTerm() = default;
    ;
};

class OneCenterTerm {
  public:
    virtual std::unique_ptr<OneCenterTerm> clone() const = 0;
    virtual void construct(
        lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a) const = 0;
    virtual ~OneCenterTerm() = default;
    ;
};

class TwoCenterTerm {
  public:
    virtual std::unique_ptr<TwoCenterTerm> clone() const = 0;
    virtual void construct(
        lexicographic::SparseMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a,
        uint32_t center_b) const = 0;
    virtual std::shared_ptr<const DenseMatrix> get_parameters() const = 0;
    virtual ~TwoCenterTerm() = default;
    ;
};
}  // namespace model::operators

#endif  //SPINNER_TERM_H
