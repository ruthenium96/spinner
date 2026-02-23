#ifndef SPINNER_EIGENFACTORIES_H
#define SPINNER_EIGENFACTORIES_H

#include "src/linalg_structures/AbstractFactories.h"

namespace spinner::linalg_structures::eigen {

class EigenDenseTransformAndDiagonalizeFactory: public AbstractDenseTransformAndDiagonalizeFactory {
  public:
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createDenseDiagonalizableMatrix(uint32_t matrix_in_space_basis_size_i) override;
    std::unique_ptr<AbstractDiagonalizableMatrix>
    createSparseDiagonalizableMatrix(uint32_t size) override;
    std::vector<std::unique_ptr<AbstractDenseVector>> createRandomUnitVectors(uint32_t size_of_vector, uint32_t number_of_vectors) override;
    std::unique_ptr<AbstractDenseVector> createVector() override;
};

} // namespace spinner::linalg_structures::eigen
#endif  //SPINNER_EIGENFACTORIES_H
