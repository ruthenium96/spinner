#ifndef SPINNER_ARMADENSEFACTORY_H
#define SPINNER_ARMADENSEFACTORY_H

#include "src/entities/data_structures/AbstractDenseFactory.h"

namespace quantum::linear_algebra {
class ArmaDenseFactory: public AbstractDenseFactory {
  public:
    std::unique_ptr<AbstractDenseMatrix> createMatrix(
        uint32_t matrix_in_space_basis_size_i,
        uint32_t matrix_in_space_basis_size_j) override;
    std::unique_ptr<AbstractDenseVector> createVector() override;
    std::vector<double>
    concatenate(const std::vector<std::unique_ptr<AbstractDenseVector>>& vectors) override;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMADENSEFACTORY_H
