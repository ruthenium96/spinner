#ifndef SPINNER_EIGENKRYLOVDENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_EIGENKRYLOVDENSESEMIUNITARYTRANSFORMER_H

#include <memory>
#include "src/linalg_structures/AbstractDenseSemiunitaryTransformer.h"
#include "EigenKrylovDenseSemiunitaryMatrix.h"

namespace spinner::linalg_structures::eigen {
template <typename T>
class EigenKrylovDenseSemiunitaryTransformer: public AbstractDenseSemiunitaryTransformer {
public:
    explicit EigenKrylovDenseSemiunitaryTransformer(const EigenKrylovDenseSemiunitaryMatrix<T>* unitary_matrix);

    std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const override;

private:
    const EigenKrylovDenseSemiunitaryMatrix<T>* unitary_matrix_;
};
    
} // namespace spinner::linalg_structures::eigen
#endif  //SPINNER_EIGENKRYLOVDENSESEMIUNITARYTRANSFORMER_H
