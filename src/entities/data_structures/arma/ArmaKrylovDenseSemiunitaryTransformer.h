#ifndef SPINNER_ARMAKRYLOVDENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_ARMAKRYLOVDENSESEMIUNITARYTRANSFORMER_H

#include <memory>
#include "src/entities/data_structures/AbstractDenseSemiunitaryTransformer.h"
#include "ArmaKrylovDenseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class ArmaKrylovDenseSemiunitaryTransformer: public AbstractDenseSemiunitaryTransformer {
public:
    explicit ArmaKrylovDenseSemiunitaryTransformer(const ArmaKrylovDenseSemiunitaryMatrix<T>* unitary_matrix);

    std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const override;

private:
    const ArmaKrylovDenseSemiunitaryMatrix<T>* unitary_matrix_;
};
    
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMAKRYLOVDENSESEMIUNITARYTRANSFORMER_H
