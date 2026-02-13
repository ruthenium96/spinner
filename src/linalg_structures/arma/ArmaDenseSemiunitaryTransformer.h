#ifndef SPINNER_ARMADENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_ARMADENSESEMIUNITARYTRANSFORMER_H

#include <memory>
#include "src/linalg_structures/AbstractDenseSemiunitaryTransformer.h"
#include "ArmaDenseSemiunitaryMatrix.h"

namespace spinner::linalg_structures {
template <typename T>
class ArmaDenseSemiunitaryTransformer: public AbstractDenseSemiunitaryTransformer {
public:
    explicit ArmaDenseSemiunitaryTransformer(const ArmaDenseSemiunitaryMatrix<T>* unitary_matrix);

    std::unique_ptr<AbstractDenseVector> calculateUnitaryTransformationOfMatrix(
        std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const override;

private:
    const ArmaDenseSemiunitaryMatrix<T>* unitary_matrix_;

    arma::Mat<T> getFirstMultiplicationResult(const AbstractDiagonalizableMatrix*) const;
    arma::Col<T> multiplyAndReturnMainDiagonal(const arma::Mat<T>& left, const arma::Mat<T>& right) const;
};
    
} // namespace spinner::linalg_structures
#endif  //SPINNER_ARMADENSESEMIUNITARYTRANSFORMER_H
