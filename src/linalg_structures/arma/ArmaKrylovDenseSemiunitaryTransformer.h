#ifndef SPINNER_ARMAKRYLOVDENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_ARMAKRYLOVDENSESEMIUNITARYTRANSFORMER_H

#include <memory>
#include "src/linalg_structures/AbstractDenseSemiunitaryTransformer.h"
#include "ArmaKrylovDenseSemiunitaryMatrix.h"

namespace spinner::linalg_structures {
template <typename T>
class ArmaKrylovDenseSemiunitaryTransformer: public AbstractDenseSemiunitaryTransformer {
public:
    explicit ArmaKrylovDenseSemiunitaryTransformer(const ArmaKrylovDenseSemiunitaryMatrix<T>* unitary_matrix);

    std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const override;

private:
    const ArmaKrylovDenseSemiunitaryMatrix<T>* unitary_matrix_;
};
    
} // namespace spinner::linalg_structures
#endif  //SPINNER_ARMAKRYLOVDENSESEMIUNITARYTRANSFORMER_H
