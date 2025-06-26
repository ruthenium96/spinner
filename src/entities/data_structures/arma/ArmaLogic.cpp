#include "ArmaLogic.h"
#include <stdexcept>

#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseDiagonalizableMatrix.h"

namespace {

template <typename T, typename M>
std::pair<arma::Col<T>, arma::Col<T>> krylovProcedureDiagonalizeValues_(
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

template <typename T, typename M>
std::tuple<arma::Col<T>, arma::Mat<T>, arma::Col<T>> krylovProcedureDiagonalizeValuesVectors_(
    const M& matrix,
    const arma::Col<T>& seed_vector,
    size_t krylov_subspace_size) {
    arma::Mat<T> krylov_matrix(krylov_subspace_size, krylov_subspace_size);
        
    if (krylov_subspace_size > matrix.n_cols) {
        throw std::invalid_argument("krylov_subspace_size bigger than size of matrix!");
    }

    arma::Mat<T> krylov_vectors(seed_vector.size(), krylov_subspace_size, arma::fill::zeros);
    krylov_vectors.col(0) = seed_vector;

    double vec_norm = arma::vecnorm(krylov_vectors.col(0));
    if (vec_norm < std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument("Extremely small value of vector norm in Krylov procedure: " + std::to_string(vec_norm));
    }
    krylov_vectors.col(0) /= vec_norm;

    auto next_vector = arma::Col<T>(seed_vector.size(), arma::fill::zeros);

    double diag_element, non_diag_element;

    for (size_t k = 0; k < krylov_subspace_size; ++k) {
        next_vector = matrix * krylov_vectors.col(k);
        diag_element = arma::dot(krylov_vectors.col(k), next_vector);
        if (k > 0) {
            non_diag_element = arma::dot(krylov_vectors.col(k-1), next_vector);
        }

        krylov_matrix.at(k, k) = diag_element;
        if (k > 0) {
            krylov_matrix.at(k-1, k) = non_diag_element;
            krylov_matrix.at(k, k-1) = non_diag_element;
        }

        next_vector -= krylov_vectors.col(k) * diag_element;
        if (k > 0) {
            next_vector -= krylov_vectors.col(k-1) * non_diag_element;
        }

        double vec_norm = arma::vecnorm(next_vector);
        if (vec_norm < std::numeric_limits<double>::epsilon()) {
            throw std::invalid_argument("Extremely small value of vector norm in Krylov procedure: " + std::to_string(vec_norm));
        }
        next_vector /= vec_norm;

        if (k + 1 < krylov_subspace_size) {
            krylov_vectors.col(k+1) = next_vector;
        }
    }

    arma::Col<T> eigenvalues;
    arma::Mat<T> eigenvectors;
    arma::eig_sym(
        eigenvalues,
        eigenvectors,
        krylov_matrix);
    
    arma::Col<T> back_projection = eigenvectors.row(0).t();

    arma::Mat<T> total_eigenvectors = krylov_vectors * eigenvectors;

    return {std::move(eigenvalues), std::move(total_eigenvectors), std::move(back_projection)};
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

template <typename T, typename M>
inline KrylovCouple krylovDiagonalizeValues_(
    const M& diagonalizableMatrix,
    const ArmaDenseVector<T>& seed_vector,
    size_t krylov_subspace_size) {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector<T>>();
    auto squared_back_projection_ = std::make_unique<ArmaDenseVector<T>>();

    auto pair = krylovProcedureDiagonalizeValues_(
        diagonalizableMatrix,
        seed_vector.getDenseVector(),
        krylov_subspace_size);

    eigenvalues_->modifyDenseVector() = std::move(pair.first);
    squared_back_projection_->modifyDenseVector() = std::move(arma::square(pair.second));

    KrylovCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.squared_back_projection = std::move(squared_back_projection_);
    return answer;
}

template <typename T>
KrylovCouple ArmaLogic<T>::krylovDiagonalizeValues(
    const AbstractDiagonalizableMatrix& diagonalizableMatrix,
    const AbstractDenseVector& seed_vector,
    size_t krylov_subspace_size) const {
    if (auto maybeDenseVector =
        dynamic_cast<const ArmaDenseVector<T>*>(&seed_vector)) {
        if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
            // It is slow, avoid this branch.
            return krylovDiagonalizeValues_(
                maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
                *maybeDenseVector,
                krylov_subspace_size
            );
        } else if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
            return krylovDiagonalizeValues_(
                maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
                *maybeDenseVector,
                krylov_subspace_size
            );
       } else {
           throw std::bad_cast();
       }
   } else {
       throw std::bad_cast();
   }
}

template <typename T, typename M>
inline KrylovTriple krylovDiagonalizeValuesVectors_(
    const M& diagonalizableMatrix,
    const ArmaDenseVector<T>& seed_vector,
    size_t krylov_subspace_size) {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector<T>>();
    auto eigenvectors_ = std::make_unique<ArmaKrylovDenseSemiunitaryMatrix<T>>();
    auto squared_back_projection_ = std::make_unique<ArmaDenseVector<T>>();

    auto triple = krylovProcedureDiagonalizeValuesVectors_(
        diagonalizableMatrix,
        seed_vector.getDenseVector(),
        krylov_subspace_size);

    eigenvalues_->modifyDenseVector() = std::move(std::get<0>(triple));
    squared_back_projection_->modifyDenseVector() = arma::square(std::get<2>(triple));
    eigenvectors_->modifyKrylovDenseSemiunitaryMatrix() = std::move(std::get<1>(triple));
    eigenvectors_->modifyBackProjectionVector() = std::move(std::get<2>(triple));
    eigenvectors_->modifySeedVector() = seed_vector.getDenseVector();

    KrylovTriple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);
    answer.squared_back_projection = std::move(squared_back_projection_);
    return answer;
}

template <typename T>
KrylovTriple ArmaLogic<T>::krylovDiagonalizeValuesVectors(
    const AbstractDiagonalizableMatrix& diagonalizableMatrix,
    const AbstractDenseVector& seed_vector,
    size_t krylov_subspace_size) const { 
    if (auto maybeDenseVector =
        dynamic_cast<const ArmaDenseVector<T>*>(&seed_vector)) {
        if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
            // It is slow, avoid this branch.
            return krylovDiagonalizeValuesVectors_(
                maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
                *maybeDenseVector,
                krylov_subspace_size
            );
        } else if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
            return krylovDiagonalizeValuesVectors_(
                maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
                *maybeDenseVector,
                krylov_subspace_size
            );
        } else {
            throw std::bad_cast();
        }
    } else {
        throw std::bad_cast();
    }
}

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
std::unique_ptr<AbstractDenseVector> ArmaLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) const {
    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return std::move(unitaryTransformAndReturnMainDiagonal_(
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
            denseSemiunitaryMatrix
        ));
    }
    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return std::move(unitaryTransformAndReturnMainDiagonal_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            denseSemiunitaryMatrix
        ));
    }
    throw std::bad_cast();
}

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
    main_diagonal->modifyDenseVector() %= denseKrylovSemiunitaryMatrix.getBackProjectionVector();
    return std::move(main_diagonal);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const ArmaKrylovDenseSemiunitaryMatrix<T>& denseKrylovSemiunitaryMatrix) const {
    if (auto maybeSparseSymmetricMatrix =
        dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(), 
            denseKrylovSemiunitaryMatrix);
    }
    if (auto maybeDenseSymmetricMatrix =
        dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(), 
            denseKrylovSemiunitaryMatrix);
    }
    throw std::bad_cast();
}

template<typename T, typename M>
inline std::unique_ptr<AbstractDiagonalizableMatrix> unitaryTransform_(
    const M& symmetric_matrix,
    const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) {
    auto result = std::make_unique<ArmaDenseDiagonalizableMatrix<T>>();

    result->modifyDenseDiagonalizableMatrix() =
        denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().t()
        * symmetric_matrix
        * denseSemiunitaryMatrix.getDenseSemiunitaryMatrix();

    return std::move(result);
}

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> ArmaLogic<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const ArmaDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) const {
    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return std::move(unitaryTransform_(
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
            denseSemiunitaryMatrix));
    }
    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return std::move(unitaryTransform_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            denseSemiunitaryMatrix));
    }
    throw std::bad_cast();
}

template class ArmaLogic<double>;
template class ArmaLogic<float>;
}  // namespace quantum::linear_algebra