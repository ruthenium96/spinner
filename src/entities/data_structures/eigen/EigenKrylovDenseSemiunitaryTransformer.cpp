#include "EigenKrylovDenseSemiunitaryTransformer.h"
#include "EigenDenseDiagonalizableMatrix.h"
#include "EigenSparseDiagonalizableMatrix.h"
#include "EigenDenseVector.h"

namespace quantum::linear_algebra {

template<typename T, typename M>
inline std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal_(
    const M& symmetric_matrix,
    const EigenKrylovDenseSemiunitaryMatrix<T>& denseKrylovSemiunitaryMatrix) {
    auto main_diagonal = std::make_unique<EigenDenseVector<T>>();

    Eigen::Vector<T, -1> firstMultiplicationResult =
    symmetric_matrix * denseKrylovSemiunitaryMatrix.getSeedVector();
    main_diagonal->modifyDenseVector() = 
        denseKrylovSemiunitaryMatrix.getKrylovDenseSemiunitaryMatrix().transpose()
        * firstMultiplicationResult;
    // instead of <n|A|r><r|n>, here we are using <n|A|r>/<r|n>,
    // putting |<r|n>|^2 in the weight of state
    main_diagonal->modifyDenseVector().array() /= denseKrylovSemiunitaryMatrix.getBackProjectionVector().array();
    return std::move(main_diagonal);
}

template <typename T>
EigenKrylovDenseSemiunitaryTransformer<T>::EigenKrylovDenseSemiunitaryTransformer(const EigenKrylovDenseSemiunitaryMatrix<T>* unitary_matrix) 
: unitary_matrix_(unitary_matrix){};

template <typename T>
std::unique_ptr<AbstractDenseVector> EigenKrylovDenseSemiunitaryTransformer<T>::calculateUnitaryTransformationOfMatrix(
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const {
    if (auto maybeSparseSymmetricMatrix =
        dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix(), 
            *unitary_matrix_);
    }
    if (auto maybeDenseSymmetricMatrix =
        dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(), 
            *unitary_matrix_);
    }
    throw std::bad_cast();    
}

template class EigenKrylovDenseSemiunitaryTransformer<double>;
template class EigenKrylovDenseSemiunitaryTransformer<float>;
} // namespace quantum::linear_algebra