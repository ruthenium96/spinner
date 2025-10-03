#ifndef SPINNER_ARMADENSESEMIUNITARYTRANSFORMER_H
#define SPINNER_ARMADENSESEMIUNITARYTRANSFORMER_H

#include <memory>
#include "src/entities/data_structures/AbstractDenseSemiunitaryTransformer.h"
#include "ArmaDenseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {
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
    
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMADENSESEMIUNITARYTRANSFORMER_H
