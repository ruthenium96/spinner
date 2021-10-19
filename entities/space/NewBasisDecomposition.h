#ifndef JULY_NEWBASISDECOMPOSITION_H
#define JULY_NEWBASISDECOMPOSITION_H

#include <cstdint>
#include <memory>
#include <ostream>

class NewBasisDecomposition {
public:
    /*
     This Iterator iterates over vectors.
     I do not know any possibilities to pass
     index_of_vector to range-base loop,
     so currently it is the best solution.
     TODO: can we pass index_of_vector to range-base loop?
     */
    struct Iterator {
        struct IndexValueItem{
            uint32_t index;
            double value;
        };
        [[nodiscard]] virtual bool hasNext() const = 0;
        virtual IndexValueItem getNext() = 0;
        virtual ~Iterator() = default;
    };

    NewBasisDecomposition();
    NewBasisDecomposition(const NewBasisDecomposition&) = delete;
    NewBasisDecomposition& operator=(const NewBasisDecomposition&) = delete;
    NewBasisDecomposition(NewBasisDecomposition&&) noexcept;
    NewBasisDecomposition& operator=(NewBasisDecomposition&&) noexcept;
    ~NewBasisDecomposition();

    [[nodiscard]] std::unique_ptr<NewBasisDecomposition::Iterator> GetNewIterator(size_t index_of_vector) const;

    /*
     Some implementations want to know the maximum size of vector.
     I pass the tensor_size to an object after all constructions, it is inconvenient.
     TODO: refactor it.
     */
    uint32_t tensor_size = 0;

    [[nodiscard]] uint32_t size() const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool vempty(uint32_t index_of_vector) const;
    void clear();

    void erase_if_zero();

    [[nodiscard]] bool is_zero(uint32_t i, uint32_t j) const;
    void move_vector_from(uint32_t i, NewBasisDecomposition& subspace_from);
    void move_all_from(NewBasisDecomposition& subspace_from);
    void copy_vector_from(uint32_t i, const NewBasisDecomposition& subspace_from);
    void copy_all_from(const NewBasisDecomposition& subspace_from);
    void resize(uint32_t new_size);

    void add_to_position(double value, uint32_t i, uint32_t j);
    double operator()(uint32_t i, uint32_t j) const;

    void normalize();

    friend std::ostream &operator<<(std::ostream &os, const NewBasisDecomposition &decomposition);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};


#endif //JULY_NEWBASISDECOMPOSITION_H
