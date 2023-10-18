#ifndef SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H
#define SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H

#include <memory>

#include "AbstractDiagonalizableMatrix.h"
#include "AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {

class AbstractSparseSemiunitaryMatrix {
  public:
    virtual ~AbstractSparseSemiunitaryMatrix() = default;
    /*
     This Iterator iterates over vectors.
     I do not know any possibilities to pass
     index_of_vector to range-base loop,
     so currently it is the best solution.
     */
    struct Iterator {
        struct IndexValueItem {
            uint32_t index;
            double value;
        };
        virtual bool hasNext() const = 0;
        virtual IndexValueItem getNext() = 0;
        virtual size_t size() const = 0;
        virtual ~Iterator() = default;
    };

    virtual std::unique_ptr<Iterator> GetNewIterator(size_t index_of_vector) const = 0;

    virtual uint32_t size_rows() const = 0;
    virtual uint32_t size_cols() const = 0;
    virtual bool empty() const = 0;
    virtual bool vempty(uint32_t index_of_vector) const = 0;
    virtual void clear() = 0;

    virtual void eraseExplicitZeros() = 0;

    virtual bool is_zero(uint32_t i, uint32_t j) const = 0;
    virtual void move_vector_from(
        uint32_t i,
        std::unique_ptr<AbstractSparseSemiunitaryMatrix>& subspace_from) = 0;

    virtual void add_to_position(double value, uint32_t i, uint32_t j) = 0;
    virtual double at(uint32_t i, uint32_t j) const = 0;

    virtual void normalize() = 0;
    virtual void unitaryTransform(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
        std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrixToAdd) const = 0;

    virtual void print(std::ostream& os) const = 0;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H
