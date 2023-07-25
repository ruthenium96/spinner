#ifndef SPINNER_ABSTRACTFACTORIES_H
#define SPINNER_ABSTRACTFACTORIES_H

#include "AbstractDenseSemiunitaryMatrix.h"
#include "AbstractSparseSemiunitaryMatrix.h"
#include "AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {
class AbstractSymmetricMatrixFactory {
  public:
    static std::shared_ptr<AbstractSymmetricMatrixFactory> defaultFactory();
    virtual std::unique_ptr<AbstractSymmetricMatrix> createDenseSymmetricMatrix(uint32_t size) = 0;
    virtual std::unique_ptr<AbstractSymmetricMatrix> createSparseSymmetricMatrix(uint32_t size) = 0;
    // TODO: move it to another class?
    virtual std::unique_ptr<AbstractDenseSemiunitaryMatrix>
    createDenseSemiunitaryMatrix(uint32_t cols, uint32_t rows) = 0;

    ~AbstractSymmetricMatrixFactory() = default;
};

class AbstractSparseSemiunitaryFactory {
  public:
    static std::shared_ptr<AbstractSparseSemiunitaryFactory> defaultSparseFactory();
    virtual std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) = 0;
};

class AbstractDenseVectorFactory {
  public:
    static std::shared_ptr<AbstractDenseVectorFactory> defaultFactory();
    virtual std::unique_ptr<AbstractDenseVector> createVector() = 0;

    ~AbstractDenseVectorFactory() = default;
};

}  // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTFACTORIES_H
