#ifndef JULY_SUBSPACE_H
#define JULY_SUBSPACE_H

#include <armadillo>
#include <cstdint>
#include <map>
#include <vector>
#include <ostream>
#include "BlockProperties.h"

#include <memory>

//#define VECTOR_OF_MAPS
// #define VECTOR_OF_SPARSE_VECTORS
//#define SPARSE_MATRIX


// #define INDEX(p) p->first
// #define VALUE(p) p->second
// using in_col_iterator = typename std::map<uint32_t, double>::iterator;
// using const_in_col_iterator = typename std::map<uint32_t, double>::const_iterator;



// TODO: create class Decomposition, move all methods to it

struct Subspace {

    struct SubspaceIterator{
        struct IndexValueItem{
            uint32_t index;
            double value;
        };
        virtual void setIndexValue(IndexValueItem) = 0;
        virtual bool hasNext() const = 0;
        virtual IndexValueItem getNext() = 0;
        virtual ~SubspaceIterator() = 0;
    };

    std::unique_ptr<SubspaceIterator> GetNewIterator();

    BlockProperties properties;
    uint32_t tensor_size = 0;

    [[nodiscard]] uint32_t size() const;
    [[nodiscard]] bool empty() const;
    [[nodiscard]] bool vempty(uint32_t index_of_vector) const;
    void clear();

    void erase_if_zero();

    [[nodiscard]] bool is_zero(uint32_t i, uint32_t j) const;
    void move_vector_from(uint32_t i, Subspace& subspace_from);
    void move_all_from(Subspace& subspace_from);
    void copy_vector_from(uint32_t i, const Subspace& subspace_from);
    void copy_all_from(const Subspace& subspace_from);
    void resize(uint32_t new_size);

    void add_to_position(double value, uint32_t i, uint32_t j);
    double operator()(uint32_t i, uint32_t j) const;


    friend std::ostream &operator<<(std::ostream &os, const Subspace &subspace);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

#endif // JULY_SUBSPACE_H
