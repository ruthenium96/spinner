#ifndef SPINNER_ABSTRACTDENSEFACTORY_H
#define SPINNER_ABSTRACTDENSEFACTORY_H

#include "AbstractDenseMatrix.h"
#include "AbstractDenseVector.h"

namespace quantum::linear_algebra {
class AbstractDenseFactory {
  public:
    static std::shared_ptr<AbstractDenseFactory> defaultFactory();
    virtual std::unique_ptr<AbstractDenseMatrix> createMatrix() = 0;
    virtual std::unique_ptr<AbstractDenseVector> createVector() = 0;
    // TODO: it is a temporary solution, fix it
    virtual std::vector<double>
    concatenate(const std::vector<std::unique_ptr<AbstractDenseVector>>& vectors) = 0;

    ~AbstractDenseFactory() = default;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTDENSEFACTORY_H
