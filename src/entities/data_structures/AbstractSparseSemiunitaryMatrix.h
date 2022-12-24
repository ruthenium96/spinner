#ifndef SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H
#define SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H

#include <memory>

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
     TODO: can we pass index_of_vector to range-base loop?
     */
    struct Iterator {
        struct IndexValueItem {
            uint32_t index;
            double value;
        };
        virtual bool hasNext() const = 0;
        virtual IndexValueItem getNext() = 0;
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
    // TODO: can also be moved to Factory
    virtual void resize(uint32_t cols, uint32_t rows) = 0;

    virtual void add_to_position(double value, uint32_t i, uint32_t j) = 0;
    virtual double at(uint32_t i, uint32_t j) const = 0;

    virtual void normalize() = 0;
    // TODO: move it from here.
    virtual void unitaryTransform(
        const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToTransform,
        std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrixToAdd) {
        size_t matrix_in_space_basis_size = this->size_cols();

        for (uint32_t index_of_space_vector_i = 0;
             index_of_space_vector_i < matrix_in_space_basis_size;
             ++index_of_space_vector_i) {
            auto outer_iterator = this->GetNewIterator(index_of_space_vector_i);
            while (outer_iterator->hasNext()) {
                auto outer_item = outer_iterator->getNext();
                uint32_t index_of_lexicographic_vector_k = outer_item.index;
                // TODO: can we start index_of_space_vector_j from index_of_space_vector_i?
                for (uint32_t index_of_space_vector_j = 0;
                     index_of_space_vector_j < matrix_in_space_basis_size;
                     ++index_of_space_vector_j) {
                    auto inner_iterator = this->GetNewIterator(index_of_space_vector_j);
                    while (inner_iterator->hasNext()) {
                        auto inner_item = inner_iterator->getNext();
                        uint32_t index_of_lexicographic_vector_l = inner_item.index;
                        double value_in_matrix_in_lexicografical_basis =
                            symmetricMatrixToTransform->at(
                                index_of_lexicographic_vector_k,
                                index_of_lexicographic_vector_l);
                        if (value_in_matrix_in_lexicografical_basis != 0) {
                            // \Delta B_{ij} = \sum_{kl} U_{ki} * \Delta A_{kl} * U_{lj}
                            symmetricMatrixToAdd->add_to_position(
                                outer_item.value * value_in_matrix_in_lexicografical_basis
                                    * inner_item.value,
                                index_of_space_vector_i,
                                index_of_space_vector_j);
                        }
                    }
                }
            }
        }
    };

    virtual void print(std::ostream& os) const = 0;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTSPARSESEMIUNITARYMATRIX_H
