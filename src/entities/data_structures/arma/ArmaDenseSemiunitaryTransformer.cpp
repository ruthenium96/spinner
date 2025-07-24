#include "ArmaDenseSemiunitaryTransformer.h"
#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaSparseDiagonalizableMatrix.h"
#include "src/entities/data_structures/arma/ArmaDenseVector.h"

namespace quantum::linear_algebra {
template <typename T, typename M>
inline std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal_(
    const M& symmetric_matrix,
    const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) {
    auto main_diagonal = std::make_unique<ArmaDenseVector<T>>();

    arma::Mat<T> firstMultiplicationResult =
        denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().t()
        * symmetric_matrix;
    main_diagonal->modifyDenseVector() = arma::diagvec(
        firstMultiplicationResult * denseSemiunitaryMatrix.getDenseSemiunitaryMatrix());
    return std::move(main_diagonal);
}

template <typename T>
ArmaDenseSemiunitaryTransformer<T>::ArmaDenseSemiunitaryTransformer(const ArmaDenseSemiunitaryMatrix<T>* unitary_matrix) 
: unitary_matrix_(unitary_matrix){};

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaDenseSemiunitaryTransformer<T>::calculateUnitaryTransformationOfMatrix(
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const {
    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return std::move(unitaryTransformAndReturnMainDiagonal_(
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
            *unitary_matrix_
        ));
    }
    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return std::move(unitaryTransformAndReturnMainDiagonal_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            *unitary_matrix_
        ));
    }
    throw std::bad_cast();
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaDenseSemiunitaryTransformer<T>::calculateUnitaryTransformationOfMatrixProduct(
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> left_matrix,
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> right_matrix) const {
    return nullptr;
}

template class ArmaDenseSemiunitaryTransformer<double>;
template class ArmaDenseSemiunitaryTransformer<float>;
} // namespace quantum::linear_algebra