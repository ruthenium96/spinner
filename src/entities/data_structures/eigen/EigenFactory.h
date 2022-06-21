#ifndef SPINNER_EIGENFACTORY_H
#define SPINNER_EIGENFACTORY_H

#include "src/entities/data_structures/AbstractFactory.h"

namespace quantum::linear_algebra {

class EigenFactory: public AbstractFactory {
  public:
    std::unique_ptr<AbstractMatrix> createMatrix() override;
    std::unique_ptr<AbstractVector> createVector() override;
    std::vector<double>
    concatenate(const std::vector<std::unique_ptr<AbstractVector>>& vectors) override;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENFACTORY_H
