#include "ArmaLogic.h"
#include <stdexcept>

#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseDiagonalizableMatrix.h"

namespace {

template <typename T, typename M>
std::pair<arma::Col<T>, arma::Col<T>> krylovDiagonalizeValuesSparse(
    const M& matrix,
    const arma::Col<T>& seed_vector,
    size_t krylov_subspace_size) {
    arma::Mat<T> krylov_matrix(krylov_subspace_size, krylov_subspace_size);
        
    if (krylov_subspace_size > matrix.n_cols) {
        throw std::invalid_argument("krylov_subspace_size bigger than size of matrix!");
    }

    auto current_vector = arma::Col<T>(seed_vector);
    double vec_norm = arma::vecnorm(current_vector);
    if (vec_norm < std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument("Extremely small value of vector norm in Krylov procedure: " + std::to_string(vec_norm));
    }
    current_vector /= vec_norm;

    auto previous_vector = arma::Col<T>(current_vector.size(), arma::fill::zeros);
    auto next_vector = arma::Col<T>(current_vector.size(), arma::fill::zeros);

    double diag_element, non_diag_element;

    for (size_t k = 0; k < krylov_subspace_size; ++k) {
        next_vector = matrix * current_vector;
        diag_element = arma::dot(current_vector, next_vector);
        if (k > 0) {
            non_diag_element = arma::dot(previous_vector, next_vector);
        }

        krylov_matrix.at(k, k) = diag_element;
        if (k > 0) {
            krylov_matrix.at(k-1, k) = non_diag_element;
            krylov_matrix.at(k, k-1) = non_diag_element;
        }

        next_vector -= current_vector * diag_element;
        if (k > 0) {
            next_vector -= previous_vector * non_diag_element;
        }

        double vec_norm = arma::vecnorm(next_vector);
        if (vec_norm < std::numeric_limits<double>::epsilon()) {
            throw std::invalid_argument("Extremely small value of vector norm in Krylov procedure: " + std::to_string(vec_norm));
        }
        next_vector /= vec_norm;

        arma::swap(previous_vector, current_vector);
        arma::swap(current_vector, next_vector);
    }

    arma::Col<T> eigenvalues;
    arma::Mat<T> eigenvectors;
    arma::eig_sym(
        eigenvalues,
        eigenvectors,
        krylov_matrix);
    
    arma::Col<T> back_projection = eigenvectors.row(0).t();

    return {std::move(eigenvalues), std::move(back_projection)};
}

} // namespace

namespace quantum::linear_algebra {

template <typename T>
std::unique_ptr<AbstractDenseVector>
ArmaLogic<T>::diagonalizeValues(const AbstractDiagonalizableMatrix& diagonalizableMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector<T>>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        auto sparseToDenseCopy = arma::Mat<T>(maybeSparseSymmetricMatrix->getSparseSymmetricMatrix());
        arma::eig_sym(eigenvalues_->modifyDenseVector(), sparseToDenseCopy);
    } else {
        throw std::bad_cast();
    }

    return eigenvalues_;
}

template <typename T>
EigenCouple ArmaLogic<T>::diagonalizeValuesVectors(
    const AbstractDiagonalizableMatrix& diagonalizableMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector<T>>();
    auto eigenvectors_ = std::make_unique<ArmaDenseSemiunitaryMatrix<T>>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
        auto sparseToDenseCopy = arma::Mat<T>(maybeSparseSymmetricMatrix->getSparseSymmetricMatrix());
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            sparseToDenseCopy);
    } else {
        throw std::bad_cast();
    }

    EigenCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);
    return answer;
}

template <typename T>
KrylovCouple ArmaLogic<T>::krylovDiagonalizeValues(
    const AbstractDiagonalizableMatrix& diagonalizableMatrix,
    const AbstractDenseVector& seed_vector,
    size_t krylov_subspace_size) const {
    
    auto eigenvalues_ = std::make_unique<ArmaDenseVector<T>>();
    auto squared_back_projection_ = std::make_unique<ArmaDenseVector<T>>();
    
    if (auto maybeDenseVector =
        dynamic_cast<const ArmaDenseVector<T>*>(&seed_vector)) {
        if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
            // It is slow, avoid this branch.
            auto pair = krylovDiagonalizeValuesSparse(
                maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
                maybeDenseVector->getDenseVector(),
                krylov_subspace_size);
            eigenvalues_->modifyDenseVector() = std::move(pair.first);
            squared_back_projection_->modifyDenseVector() = std::move(arma::square(pair.second));
        } else if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {

            auto pair = krylovDiagonalizeValuesSparse(
                maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
                maybeDenseVector->getDenseVector(),
                krylov_subspace_size);
            eigenvalues_->modifyDenseVector() = std::move(pair.first);
            squared_back_projection_->modifyDenseVector() = std::move(arma::square(pair.second));
       } else {
           throw std::bad_cast();
       }
   } else {
       throw std::bad_cast();
   }

    KrylovCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.squared_back_projection = std::move(squared_back_projection_);
    return answer;
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const std::unique_ptr<AbstractDenseSemiunitaryMatrix>& denseSemiunitaryMatrix) const {
    return unitaryTransformAndReturnMainDiagonal(symmetricMatrix, *denseSemiunitaryMatrix);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) const {
    auto main_diagonal = std::make_unique<ArmaDenseVector<T>>();

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        arma::Mat<T> firstMultiplicationResult =
            denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix();
        main_diagonal->modifyDenseVector() = arma::diagvec(
            firstMultiplicationResult * denseSemiunitaryMatrix.getDenseSemiunitaryMatrix());
        return std::move(main_diagonal);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        arma::Mat<T> firstMultiplicationResult =
            denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix();
        main_diagonal->modifyDenseVector() = arma::diagvec(
            firstMultiplicationResult * denseSemiunitaryMatrix.getDenseSemiunitaryMatrix());
        return std::move(main_diagonal);
    }

    throw std::bad_cast();
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> ArmaLogic<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) const {
    auto result = std::make_unique<ArmaDenseDiagonalizableMatrix<T>>();


    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix()
            * denseSemiunitaryMatrix.getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix()
            * denseSemiunitaryMatrix.getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    throw std::bad_cast();
}

template class ArmaLogic<double>;
template class ArmaLogic<float>;
}  // namespace quantum::linear_algebra