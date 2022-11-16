#ifndef SPINNER_ABSTRACTSPARSEMATRIX_H
#define SPINNER_ABSTRACTSPARSEMATRIX_H

#include <memory>
namespace quantum::linear_algebra {

class AbstractSparseMatrix {
  public:
    virtual ~AbstractSparseMatrix() = default;
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

    static std::unique_ptr<AbstractSparseMatrix> defaultSparseMatrix();

    virtual std::unique_ptr<Iterator> GetNewIterator(size_t index_of_vector) const = 0;

    virtual uint32_t size() const = 0;
    virtual bool empty() const = 0;
    virtual bool vempty(uint32_t index_of_vector) const = 0;
    virtual void clear() = 0;

    virtual void erase_if_zero() = 0;

    virtual bool is_zero(uint32_t i, uint32_t j) const = 0;
    virtual void
    move_vector_from(uint32_t i, std::unique_ptr<AbstractSparseMatrix>& subspace_from) = 0;
    virtual void move_all_from(std::unique_ptr<AbstractSparseMatrix>& subspace_from) = 0;
    virtual void
    copy_vector_from(uint32_t i, const std::unique_ptr<AbstractSparseMatrix>& subspace_from) = 0;
    virtual void copy_all_from(const std::unique_ptr<AbstractSparseMatrix>& subspace_from) = 0;
    virtual void resize(uint32_t new_size) = 0;

    virtual void add_to_position(double value, uint32_t i, uint32_t j) = 0;
    virtual double at(uint32_t i, uint32_t j) const = 0;

    virtual void normalize() = 0;

    virtual bool
    is_equal_up_to_vector_order(const std::unique_ptr<AbstractSparseMatrix>& rhs) const = 0;

    virtual void print(std::ostream& os) const = 0;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ABSTRACTSPARSEMATRIX_H
