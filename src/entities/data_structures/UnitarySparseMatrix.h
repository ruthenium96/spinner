#ifndef SPINNER_UNITARYSPARSEMATRIX_H
#define SPINNER_UNITARYSPARSEMATRIX_H

#include <cstdint>
#include <memory>
#include <ostream>

class UnitarySparseMatrix {
  public:
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

    UnitarySparseMatrix();
    UnitarySparseMatrix(const UnitarySparseMatrix&) = delete;
    UnitarySparseMatrix& operator=(const UnitarySparseMatrix&) = delete;
    UnitarySparseMatrix(UnitarySparseMatrix&&) noexcept;
    UnitarySparseMatrix& operator=(UnitarySparseMatrix&&) noexcept;
    ~UnitarySparseMatrix();

    std::unique_ptr<UnitarySparseMatrix::Iterator> GetNewIterator(size_t index_of_vector) const;

    /*
     Some implementations want to know the maximum size of vector.
     I pass the tensor_size to an object after all constructions, it is inconvenient.
     TODO: refactor it.
     */
    uint32_t tensor_size = 0;

    uint32_t size() const;
    bool empty() const;
    bool vempty(uint32_t index_of_vector) const;
    void clear();

    void erase_if_zero();

    bool is_zero(uint32_t i, uint32_t j) const;
    void move_vector_from(uint32_t i, UnitarySparseMatrix& subspace_from);
    void move_all_from(UnitarySparseMatrix& subspace_from);
    void copy_vector_from(uint32_t i, const UnitarySparseMatrix& subspace_from);
    void copy_all_from(const UnitarySparseMatrix& subspace_from);
    void resize(uint32_t new_size);

    void add_to_position(double value, uint32_t i, uint32_t j);
    double operator()(uint32_t i, uint32_t j) const;

    void normalize();

    bool is_equal_up_to_vector_order(const UnitarySparseMatrix& rhs) const;

    friend std::ostream& operator<<(std::ostream& os, const UnitarySparseMatrix& decomposition);

  private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

#endif  //SPINNER_UNITARYSPARSEMATRIX_H
