#include "ArmaLogic.h"

#include "ArmaDenseSemiunitaryMatrix.h"
#include "ArmaDenseSymmetricMatrix.h"
#include "ArmaDenseVector.h"
#include "ArmaSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractDenseVector> ArmaLogic::diagonalizeValues(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix) const {
    return diagonalizeValues(*symmetricMatrix);
}

std::unique_ptr<AbstractDenseVector>
ArmaLogic::diagonalizeValues(const AbstractSymmetricMatrix& symmetricMatrix) const {
    // TODO: refactor this code
    auto eigenvalues = std::make_unique<ArmaDenseVector>();

    auto maybeDenseSymmetricMatrix =
        dynamic_cast<const ArmaDenseSymmetricMatrix*>(&symmetricMatrix);
    if (maybeDenseSymmetricMatrix != nullptr) {
        arma::eig_sym(
            eigenvalues->modifyDenseVector(),
            maybeDenseSymmetricMatrix->getDenseSymmetricMatrix());
        return eigenvalues;
    }

    auto maybeSparseSymmetricMatrix =
        dynamic_cast<const ArmaSparseSymmetricMatrix*>(&symmetricMatrix);
    if (maybeSparseSymmetricMatrix != nullptr) {
        size_t size = maybeSparseSymmetricMatrix->size();
        arma::eigs_sym(
            eigenvalues->modifyDenseVector(),
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
            size - 1,
            "la");
        return eigenvalues;
    }
    throw std::bad_cast();
}

EigenCouple ArmaLogic::diagonalizeValuesVectors(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix) const {
    return diagonalizeValuesVectors(*symmetricMatrix);
}

EigenCouple
ArmaLogic::diagonalizeValuesVectors(const AbstractSymmetricMatrix& symmetricMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector>();
    auto eigenvectors_ = std::make_unique<ArmaDenseSemiunitaryMatrix>();

    auto maybeDenseSymmetricMatrix =
        dynamic_cast<const ArmaDenseSymmetricMatrix*>(&symmetricMatrix);
    if (maybeDenseSymmetricMatrix != nullptr) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            maybeDenseSymmetricMatrix->getDenseSymmetricMatrix());

        EigenCouple answer;
        answer.eigenvalues = std::move(eigenvalues_);
        answer.eigenvectors = std::move(eigenvectors_);
        return answer;
    }

    auto maybeSparseSymmetricMatrix =
        dynamic_cast<const ArmaSparseSymmetricMatrix*>(&symmetricMatrix);
    if (maybeSparseSymmetricMatrix != nullptr) {
        size_t size = maybeSparseSymmetricMatrix->size();
        arma::eigs_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            maybeSparseSymmetricMatrix->getSparseSymmetricMatrix(),
            size - 1,
            "la");

        EigenCouple answer;
        answer.eigenvalues = std::move(eigenvalues_);
        answer.eigenvectors = std::move(eigenvectors_);
        return answer;
    }

    throw std::bad_cast();
}

std::unique_ptr<AbstractDenseVector> ArmaLogic::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix,
    const std::unique_ptr<AbstractDenseSemiunitaryMatrix>& denseSemiunitaryMatrix) const {
    return unitaryTransformAndReturnMainDiagonal(symmetricMatrix, *denseSemiunitaryMatrix);
}

std::unique_ptr<AbstractDenseVector> ArmaLogic::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto main_diagonal = std::make_unique<ArmaDenseVector>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const ArmaDenseSemiunitaryMatrix*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    auto maybeSparseSymmetricMatrix =
        dynamic_cast<const ArmaSparseSymmetricMatrix*>(symmetricMatrix.get());
    if (maybeSparseSymmetricMatrix != nullptr) {
        arma::dmat transformed_matrix = maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix().t()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        main_diagonal->modifyDenseVector() = transformed_matrix.diag();
        return std::move(main_diagonal);
    }

    auto maybeDenseSymmetricMatrix =
        dynamic_cast<const ArmaDenseSymmetricMatrix*>(symmetricMatrix.get());
    if (maybeDenseSymmetricMatrix != nullptr) {
        arma::dmat transformed_matrix = maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseSymmetricMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        main_diagonal->modifyDenseVector() = transformed_matrix.diag();
        return std::move(main_diagonal);
    }

    throw std::bad_cast();
}
}  // namespace quantum::linear_algebra