#include "ArmaKrylovDenseSemiunitaryTransformer.h"
#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaSparseDiagonalizableMatrix.h"
#include "ArmaDenseVector.h"

namespace quantum::linear_algebra {

template<typename T, typename M>
inline std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal_(
    const M& symmetric_matrix,
    const ArmaKrylovDenseSemiunitaryMatrix<T>& denseKrylovSemiunitaryMatrix) {
    auto main_diagonal = std::make_unique<ArmaDenseVector<T>>();

    arma::Col<T> firstMultiplicationResult =
    symmetric_matrix * denseKrylovSemiunitaryMatrix.getSeedVector();
    main_diagonal->modifyDenseVector() = 
        denseKrylovSemiunitaryMatrix.getKrylovDenseSemiunitaryMatrix().t()
        * firstMultiplicationResult;
    // instead of <n|A|r><r|n>, here we are using <n|A|r>/<r|n>, 
    // putting |<r|n>|^2 in the weight of state
    main_diagonal->modifyDenseVector() /= denseKrylovSemiunitaryMatrix.getBackProjectionVector();
    return std::move(main_diagonal);
}

template <typename T>
ArmaKrylovDenseSemiunitaryTransformer<T>::ArmaKrylovDenseSemiunitaryTransformer(const ArmaKrylovDenseSemiunitaryMatrix<T>* unitary_matrix) 
: unitary_matrix_(unitary_matrix){};

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaKrylovDenseSemiunitaryTransformer<T>::calculateUnitaryTransformationOfMatrix(
    std::reference_wrapper<const std::unique_ptr<AbstractDiagonalizableMatrix>> matrix) const {
    if (auto maybeSparseSymmetricMatrix =
        dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(), 
            *unitary_matrix_);
    }
    if (auto maybeDenseSymmetricMatrix =
        dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(matrix.get().get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(), 
            *unitary_matrix_);
    }
    throw std::bad_cast();
}

template class ArmaKrylovDenseSemiunitaryTransformer<double>;
template class ArmaKrylovDenseSemiunitaryTransformer<float>;
} // namespace quantum::linear_algebra