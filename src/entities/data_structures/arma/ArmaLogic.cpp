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
    auto eigenvalues_ = std::make_unique<ArmaDenseVector>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseSymmetricMatrix*>(&symmetricMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            maybeDenseSymmetricMatrix->getDenseSymmetricMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseSymmetricMatrix*>(&symmetricMatrix)) {
        auto sparseToDenseCopy = arma::dmat(maybeSparseSymmetricMatrix->getSparseSymmetricMatrix());
        arma::eig_sym(eigenvalues_->modifyDenseVector(), sparseToDenseCopy);
    } else {
        throw std::bad_cast();
    }

    return eigenvalues_;
}

EigenCouple ArmaLogic::diagonalizeValuesVectors(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix) const {
    return diagonalizeValuesVectors(*symmetricMatrix);
}

EigenCouple
ArmaLogic::diagonalizeValuesVectors(const AbstractSymmetricMatrix& symmetricMatrix) const {
    auto eigenvalues_ = std::make_unique<ArmaDenseVector>();
    auto eigenvectors_ = std::make_unique<ArmaDenseSemiunitaryMatrix>();

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseSymmetricMatrix*>(&symmetricMatrix)) {
        arma::eig_sym(
            eigenvalues_->modifyDenseVector(),
            eigenvectors_->modifyDenseSemiunitaryMatrix(),
            maybeDenseSymmetricMatrix->getDenseSymmetricMatrix());
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseSymmetricMatrix*>(&symmetricMatrix)) {
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

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const ArmaSparseSymmetricMatrix*>(symmetricMatrix.get())) {
        arma::dmat firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix();
        main_diagonal->resize(maybeSparseSymmetricMatrix->size());
#pragma omp parallel for shared( \
    main_diagonal, \
    firstMultiplicationResult, \
    maybeDenseSemiunitaryMatrix) default(none)
        for (size_t i = 0; i < main_diagonal->size(); ++i) {
            auto value = firstMultiplicationResult.row(i)
                * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().col(i);
            main_diagonal->modifyDenseVector().at(i) = as_scalar(value);
        }
        return std::move(main_diagonal);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const ArmaDenseSymmetricMatrix*>(symmetricMatrix.get())) {
        arma::dmat firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().t()
            * maybeDenseSymmetricMatrix->getDenseSymmetricMatrix();
        main_diagonal->resize(maybeDenseSymmetricMatrix->size());
#pragma omp parallel for shared( \
    main_diagonal, \
    firstMultiplicationResult, \
    maybeDenseSemiunitaryMatrix) default(none)
        for (size_t i = 0; i < main_diagonal->size(); ++i) {
            auto value = firstMultiplicationResult.row(i)
                * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().col(i);
            main_diagonal->modifyDenseVector().at(i) = as_scalar(value);
        }
        return std::move(main_diagonal);
    }

    throw std::bad_cast();
}
}  // namespace quantum::linear_algebra