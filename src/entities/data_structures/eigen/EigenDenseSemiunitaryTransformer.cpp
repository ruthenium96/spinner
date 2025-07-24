#include "EigenDenseSemiunitaryTransformer.h"
#include "EigenDenseDiagonalizableMatrix.h"
#include "EigenSparseDiagonalizableMatrix.h"
#include "EigenDenseVector.h"

namespace quantum::linear_algebra {

template <typename T, typename M>
inline std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal_(
    const M& symmetricMatrix,
    const EigenDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) {
    auto main_diagonal = std::make_unique<EigenDenseVector<T>>();

    Eigen::Matrix<T, -1, -1> firstMultiplicationResult =
        denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().transpose()
        * symmetricMatrix;
    main_diagonal->resize(symmetricMatrix.rows());
#pragma omp parallel for shared( \
    main_diagonal, \
    firstMultiplicationResult, \
    denseSemiunitaryMatrix) default(none)
    for (size_t i = 0; i < main_diagonal->size(); ++i) {
        auto value = firstMultiplicationResult.row(i)
            * denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().col(i);
        main_diagonal->modifyDenseVector().coeffRef(i) = value;
    }
    return std::move(main_diagonal);
}

template <typename T>
EigenDenseSemiunitaryTransformer<T>::EigenDenseSemiunitaryTransformer(const EigenDenseSemiunitaryMatrix<T>* unitary_matrix) 
: unitary_matrix_(unitary_matrix){};

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenDenseSemiunitaryTransformer<T>::calculateUnitaryTransformationOfMatrix(
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const {
    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix().transpose(),
            *unitary_matrix_
        );
    }
    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            *unitary_matrix_
        );
    }
    throw std::bad_cast();
}

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenDenseSemiunitaryTransformer<T>::calculateUnitaryTransformationOfMatrixProduct(
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> left_matrix,
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> right_matrix) const {
    return nullptr;
}

template class EigenDenseSemiunitaryTransformer<double>;
template class EigenDenseSemiunitaryTransformer<float>;
} // namespace quantum::linear_algebra