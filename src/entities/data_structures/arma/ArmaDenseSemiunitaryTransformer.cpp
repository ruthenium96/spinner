#include "ArmaDenseSemiunitaryTransformer.h"
#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaSparseDiagonalizableMatrix.h"
#include "src/entities/data_structures/arma/ArmaDenseVector.h"

namespace quantum::linear_algebra {
template <typename T, typename M>
inline arma::Mat<T> first_multiplication_(
    const M& symmetric_matrix,
    const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) {
    return denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().t() * symmetric_matrix;
}

template <typename T>
ArmaDenseSemiunitaryTransformer<T>::ArmaDenseSemiunitaryTransformer(const ArmaDenseSemiunitaryMatrix<T>* unitary_matrix) 
: unitary_matrix_(unitary_matrix){};

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaDenseSemiunitaryTransformer<T>::calculateUnitaryTransformationOfMatrix(
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const {
    const auto& first_multiplication_result = getFirstMultiplicationResult(matrix.get().get());
    auto main_diagonal = std::make_unique<ArmaDenseVector<T>>();
    main_diagonal->modifyDenseVector() = std::move(multiplyAndReturnMainDiagonal(first_multiplication_result, unitary_matrix_->getDenseSemiunitaryMatrix()));
    return std::move(main_diagonal);
}

template <typename T>
arma::Mat<T> ArmaDenseSemiunitaryTransformer<T>::getFirstMultiplicationResult(const AbstractDiagonalizableMatrix* ptr) const {
    if (auto maybeSparseSymmetricMatrix =
        dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(ptr)) {
        return std::move(first_multiplication_(
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
            *unitary_matrix_
        ));
    } else if (auto maybeDenseSymmetricMatrix =
        dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(ptr)) {   
        return std::move(first_multiplication_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            *unitary_matrix_
        ));
    } else {
        throw std::bad_cast();
    }
}

template <typename T>
arma::Col<T> ArmaDenseSemiunitaryTransformer<T>::multiplyAndReturnMainDiagonal(const arma::Mat<T>& left, const arma::Mat<T>& right) const {
    return arma::diagvec(left * right);
}

template class ArmaDenseSemiunitaryTransformer<double>;
template class ArmaDenseSemiunitaryTransformer<float>;
} // namespace quantum::linear_algebra