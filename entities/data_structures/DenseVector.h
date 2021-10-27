#ifndef JULY_DENSEVECTOR_H
#define JULY_DENSEVECTOR_H

#include "entities/data_structures/DenseMatrix.h"
#include <memory>

class DenseVector {
    friend DenseMatrix;
public:
    DenseVector();
    DenseVector(const DenseVector&) = delete;
    DenseVector& operator=(const DenseVector&) = delete;
    DenseVector(DenseVector&&) noexcept;
    DenseVector& operator=(DenseVector&&) noexcept;
    ~DenseVector();

    friend std::ostream &operator<<(std::ostream &os, const DenseVector &raw_data);
private:
    class SubspectrumDataImpl;
    std::unique_ptr<SubspectrumDataImpl> pImpl;
};

#endif //JULY_DENSEVECTOR_H
