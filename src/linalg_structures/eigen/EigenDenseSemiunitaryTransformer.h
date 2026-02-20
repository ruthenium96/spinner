#ifndef SPINNER_EIGENDENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_EIGENDENSESEMIUNITARYTRANSFORMER_H

#include <memory>
#include "src/linalg_structures/AbstractDenseSemiunitaryTransformer.h"
#include "EigenDenseSemiunitaryMatrix.h"

namespace spinner::linalg_structures::eigen {
template <typename T>
class EigenDenseSemiunitaryTransformer: public AbstractDenseSemiunitaryTransformer {
public:
    explicit EigenDenseSemiunitaryTransformer(const EigenDenseSemiunitaryMatrix<T>* unitary_matrix);

    std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const override;

private:
    const EigenDenseSemiunitaryMatrix<T>* unitary_matrix_;
};
    
} // namespace spinner::linalg_structures::eigen
#endif  //SPINNER_EIGENDENSESEMIUNITARYTRANSFORMER_H
