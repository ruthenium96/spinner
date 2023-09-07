#include "ArmaLogic.h"

#include "ArmaDenseDiagonalizableMatrix.h"
#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseDiagonalizableMatrix.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractDenseVector> ArmaLogic::diagonalizeValues(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& diagonalizableMatrix) const {
    return diagonalizeValues(*diagonalizableMatrix);
}

std::unique_ptr<AbstractDenseVector>
ArmaLogic::diagonalizeValues(const AbstractDiagonalizableMatrix& diagonalizableMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix*>(&diagonalizableMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix*>(&diagonalizableMatrix)) {
        auto sparseToDenseCopy = arma::dmat(maybeSparseSymmetricMatrix->getSparseSymmetricMatrix());
        arma::eig_sym(eigenvalues_->modifyDenseVector(), sparseToDenseCopy);
    } else {
        throw std::bad_cast();
    }

    return eigenvalues_;
}

EigenCouple ArmaLogic::diagonalizeValuesVectors(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix) const {
    return diagonalizeValuesVectors(*symmetricMatrix);
}

EigenCouple ArmaLogic::diagonalizeValuesVectors(
    const AbstractDiagonalizableMatrix& diagonalizableMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector>();
    auto eigenvectors_ = std::make_unique<ArmaDenseSemiunitaryMatrix>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix*>(&diagonalizableMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix*>(&diagonalizableMatrix)) {
        auto sparseToDenseCopy = arma::dmat(maybeSparseSymmetricMatrix->getSparseSymmetricMatrix());
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

std::unique_ptr<AbstractDenseVector> ArmaLogic::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const std::unique_ptr<AbstractDenseSemiunitaryMatrix>& denseSemiunitaryMatrix) const {
    return unitaryTransformAndReturnMainDiagonal(symmetricMatrix, *denseSemiunitaryMatrix);
}

std::unique_ptr<AbstractDenseVector> ArmaLogic::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto main_diagonal = std::make_unique<ArmaDenseVector>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const ArmaDenseSemiunitaryMatrix*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix*>(symmetricMatrix.get())) {
        arma::dmat firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix();
        main_diagonal->modifyDenseVector() = arma::diagvec(
            firstMultiplicationResult * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix());
        return std::move(main_diagonal);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix*>(symmetricMatrix.get())) {
        arma::dmat firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix();
        main_diagonal->modifyDenseVector() = arma::diagvec(
            firstMultiplicationResult * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix());
        return std::move(main_diagonal);
    }

    throw std::bad_cast();
}

std::unique_ptr<AbstractDiagonalizableMatrix> ArmaLogic::unitaryTransform(
    const std::unique_ptr<AbstractDiagonalizableMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto result = std::make_unique<ArmaDenseDiagonalizableMatrix>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const ArmaDenseSemiunitaryMatrix*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseDiagonalizableMatrix*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseDiagonalizableMatrix*>(symmetricMatrix.get())) {
        result->modifyDenseDiagonalizableMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseDiagonalizableMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    throw std::bad_cast();
}
}  // namespace quantum::linear_algebra