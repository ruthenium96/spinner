#include "EigenLogic.h"
#include <stdexcept>

#include "EigenDenseDiagonalizableMatrix.h"
#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseVector.h"
#include "EigenSparseDiagonalizableMatrix.h"

namespace {

template <typename T, typename M>
std::pair<Eigen::Vector<T, -1>, Eigen::Vector<T, -1>> krylovProcedureDiagonalizeValues_(
    const M& matrix,
    const Eigen::Vector<T, -1>& seed_vector,
    size_t krylov_subspace_size) {
    Eigen::Matrix<T, -1, -1> krylov_matrix(krylov_subspace_size, krylov_subspace_size);
    krylov_matrix.fill(0);
    
    if (krylov_subspace_size > matrix.cols()) {
        throw std::invalid_argument("krylov_subspace_size bigger than size of matrix!");
    }

    // use objects instead of pointers and use swap?
    auto current_vector = std::make_unique<Eigen::Vector<T, -1>>(seed_vector);
    double vec_norm = current_vector->norm();
    if (vec_norm < std::numeric_limits<double>::epsilon()) {
        throw std::invalid_argument("Extremely small value of vector norm in Krylov procedure: " + std::to_string(vec_norm));
    }
    *current_vector /= vec_norm;

    auto previous_vector = std::make_unique<Eigen::Vector<T, -1>>(current_vector->size());
    previous_vector->fill(0);
    auto next_vector = std::make_unique<Eigen::Vector<T, -1>>(current_vector->size());
    next_vector->fill(0);

    double diag_element, non_diag_element;

    for (size_t k = 0; k < krylov_subspace_size; ++k) {
        *next_vector = (matrix) * (*current_vector);
        diag_element = current_vector->dot(*next_vector);
        if (k > 0) {
            non_diag_element = previous_vector->dot(*next_vector);
        }

        krylov_matrix(k, k) = diag_element;
        if (k > 0) {
            krylov_matrix(k-1, k) = non_diag_element;
            krylov_matrix(k, k-1) = non_diag_element;
        }

        *next_vector -= *current_vector * diag_element;
        if (k > 0) {
            *next_vector -= *previous_vector * non_diag_element;
        }

        double vec_norm = next_vector->norm();
        if (vec_norm < std::numeric_limits<double>::epsilon()) {
            throw std::invalid_argument("Extremely small value of vector norm in Krylov procedure: " + std::to_string(vec_norm));
        }
        *next_vector /= vec_norm;

        // use objects instead of pointers and use arma::swap?
        std::swap(previous_vector, current_vector);
        std::swap(current_vector, next_vector);
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, -1, -1>> es;
    es.compute(
        krylov_matrix,
        Eigen::ComputeEigenvectors);
    auto eigenvalues = std::move(es.eigenvalues());
    auto eigenvectors = std::move(es.eigenvectors());
        
    auto back_projection = eigenvectors.row(0);

    return {std::move(eigenvalues), std::move(back_projection)};
}   

} // namespace

namespace quantum::linear_algebra {

template <typename T, typename M>
inline std::unique_ptr<AbstractDenseVector> unitaryTransformAndReturnMainDiagonal_(
    const M& symmetricMatrix,
    const EigenDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) {
    auto main_diagonal = std::make_unique<EigenDenseVector<T>>();

    Eigen::Matrix<T, -1, -1> firstMultiplicationResult =
        denseSemiunitaryMatrix.getDenseSemiunitaryMatrix().transpose()
        * symmetricMatrix;
    main_diagonal->resize(symmetricMatrix.size());
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
std::unique_ptr<AbstractDenseVector> EigenLogic<T>::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const EigenDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) const {
    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const EigenDenseSemiunitaryMatrix<T>*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }
    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix().transpose(),
            *maybeDenseSemiunitaryMatrix
        );
    }
    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        return unitaryTransformAndReturnMainDiagonal_(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            *maybeDenseSemiunitaryMatrix
        );
    }
    throw std::bad_cast();
}

template <typename T, typename M>
inline KrylovCouple krylovDiagonalizeValues_(
    const M& diagonalizableMatrix,
    const EigenDenseVector<T>& seed_vector,
    size_t krylov_subspace_size) {
    auto eigenvalues_ = std::make_unique<EigenDenseVector<T>>();
    auto squared_back_projection_ = std::make_unique<EigenDenseVector<T>>();
    auto pair = krylovProcedureDiagonalizeValues_(
        diagonalizableMatrix,
        seed_vector.getDenseVector(),
        krylov_subspace_size);
    eigenvalues_->modifyDenseVector() = std::move(pair.first);
    squared_back_projection_->modifyDenseVector() = std::move(pair.second.array().square());

    KrylovCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.squared_back_projection = std::move(squared_back_projection_);
    return answer;
}

template <typename T>
KrylovCouple EigenLogic<T>::krylovDiagonalizeValues(
    const AbstractDiagonalizableMatrix& diagonalizableMatrix,
    const AbstractDenseVector& seed_vector,
    size_t krylov_subspace_size) const {    
    if (auto maybeDenseVector =
        dynamic_cast<const EigenDenseVector<T>*>(&seed_vector)) {
        if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
            return krylovDiagonalizeValues_(
                maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
                *maybeDenseVector,
                krylov_subspace_size
            );
        } else if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(&diagonalizableMatrix)) {
            return krylovDiagonalizeValues_(
                maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix(),
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

template <typename T>
std::unique_ptr<AbstractDiagonalizableMatrix> EigenLogic<T>::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const EigenDenseSemiunitaryMatrix<T>& denseSemiunitaryMatrix) const {
    auto result = std::make_unique<EigenDenseDiagonalizableMatrix<T>>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const EigenDenseSemiunitaryMatrix<T>*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        // TODO: check if all transpositions are efficient
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix().transpose()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    throw std::bad_cast();
}

template <typename T>
std::unique_ptr<AbstractDenseVector>
EigenLogic<T>::diagonalizeValues(const AbstractDiagonalizableMatrix& symmetricMatrix) const {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, -1, -1>> es;

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        es.compute(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            Eigen::EigenvaluesOnly);
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        Eigen::Matrix<T, -1, -1> sparseToDenseCopy =
            maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix();
        es.compute(sparseToDenseCopy, Eigen::EigenvaluesOnly);
    } else {
        throw std::bad_cast();
    }

    auto eigenvalues_ = std::make_unique<EigenDenseVector<T>>();
    eigenvalues_->modifyDenseVector() = es.eigenvalues();

    return eigenvalues_;
}

template <typename T>
EigenCouple
EigenLogic<T>::diagonalizeValuesVectors(const AbstractDiagonalizableMatrix& symmetricMatrix) const {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, -1, -1>> es;

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        es.compute(
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix(),
            Eigen::ComputeEigenvectors);
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseDiagonalizableMatrix<T>*>(&symmetricMatrix)) {
        Eigen::Matrix<T, -1, -1> sparseToDenseCopy =
            maybeSparseSymmetricMatrix->getSparseDiagonalizableMatrix();
        es.compute(sparseToDenseCopy, Eigen::ComputeEigenvectors);
    } else {
        throw std::bad_cast();
    }

    auto eigenvalues_ = std::make_unique<EigenDenseVector<T>>();
    auto eigenvectors_ = std::make_unique<EigenDenseSemiunitaryMatrix<T>>();
    eigenvalues_->modifyDenseVector() = es.eigenvalues();
    eigenvectors_->modifyDenseSemiunitaryMatrix() = es.eigenvectors();

    EigenCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);
    return answer;
}

template class EigenLogic<double>;
template class EigenLogic<float>;
}  // namespace quantum::linear_algebra