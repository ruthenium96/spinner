#ifndef SPINNER_ABSTRACTFACTORY_H
#define SPINNER_ABSTRACTFACTORY_H

#include "AbstractMatrix.h"
#include "AbstractVector.h"

namespace quantum::linear_algebra {
class AbstractFactory {
  public:
    virtual std::unique_ptr<AbstractMatrix> createMatrix() = 0;
    virtual std::unique_ptr<AbstractVector> createVector() = 0;
    // TODO: it is a temporary solution, fix it
    virtual std::vector<double>
    concatenate(const std::vector<std::unique_ptr<AbstractVector>>& vectors) = 0;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTFACTORY_H
