#ifndef SPINNER_EIGENKRYLOVDENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_EIGENKRYLOVDENSESEMIUNITARYTRANSFORMER_H

#include <memory>
#include "src/entities/data_structures/AbstractDenseSemiunitaryTransformer.h"
#include "EigenKrylovDenseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
template <typename T>
class EigenKrylovDenseSemiunitaryTransformer: public AbstractDenseSemiunitaryTransformer {
public:
    explicit EigenKrylovDenseSemiunitaryTransformer(const EigenKrylovDenseSemiunitaryMatrix<T>* unitary_matrix);

    std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const override;

private:
    const EigenKrylovDenseSemiunitaryMatrix<T>* unitary_matrix_;
};
    
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENKRYLOVDENSESEMIUNITARYTRANSFORMER_H
