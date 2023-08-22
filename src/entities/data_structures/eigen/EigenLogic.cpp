#include "EigenLogic.h"

#include "EigenDenseSemiunitaryMatrix.h"
#include "EigenDenseSymmetricMatrix.h"
#include "EigenDenseVector.h"
#include "EigenSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {

std::unique_ptr<AbstractDenseVector> EigenLogic::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto main_diagonal = std::make_unique<EigenDenseVector>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const EigenDenseSemiunitaryMatrix*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseSymmetricMatrix*>(symmetricMatrix.get())) {
        Eigen::MatrixXd firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix().transpose();
        main_diagonal->resize(maybeSparseSymmetricMatrix->size());
#pragma omp parallel for shared( \
    main_diagonal, \
    firstMultiplicationResult, \
    maybeDenseSemiunitaryMatrix) default(none)
        for (size_t i = 0; i < main_diagonal->size(); ++i) {
            auto value = firstMultiplicationResult.row(i)
                * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().col(i);
            main_diagonal->modifyDenseVector().coeffRef(i) = value;
        }
        return std::move(main_diagonal);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseSymmetricMatrix*>(symmetricMatrix.get())) {
        Eigen::MatrixXd firstMultiplicationResult =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeDenseSymmetricMatrix->getDenseSymmetricMatrix();
        main_diagonal->resize(maybeDenseSymmetricMatrix->size());
#pragma omp parallel for shared( \
        main_diagonal, \
            firstMultiplicationResult, \
            maybeDenseSemiunitaryMatrix) default(none)
        for (size_t i = 0; i < main_diagonal->size(); ++i) {
            auto value = firstMultiplicationResult.row(i)
                * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().col(i);
            main_diagonal->modifyDenseVector().coeffRef(i) = value;
        }
        return std::move(main_diagonal);
    }

    throw std::bad_cast();
}

std::unique_ptr<AbstractSymmetricMatrix> EigenLogic::unitaryTransform(
    const std::unique_ptr<AbstractSymmetricMatrix>& symmetricMatrix,
    const AbstractDenseSemiunitaryMatrix& denseSemiunitaryMatrix) const {
    auto result = std::make_unique<EigenDenseSymmetricMatrix>();

    auto maybeDenseSemiunitaryMatrix =
        dynamic_cast<const EigenDenseSemiunitaryMatrix*>(&denseSemiunitaryMatrix);
    if (maybeDenseSemiunitaryMatrix == nullptr) {
        throw std::bad_cast();
    }

    if (auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseSymmetricMatrix*>(symmetricMatrix.get())) {
        // TODO: check if all transpositions are efficient
        result->modifyDenseSymmetricMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeSparseSymmetricMatrix->getSparseSymmetricMatrix().transpose()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseSymmetricMatrix*>(symmetricMatrix.get())) {
        result->modifyDenseSymmetricMatrix() =
            maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix().transpose()
            * maybeDenseSymmetricMatrix->getDenseSymmetricMatrix()
            * maybeDenseSemiunitaryMatrix->getDenseSemiunitaryMatrix();
        return std::move(result);
    }

    throw std::bad_cast();
}

std::unique_ptr<AbstractDenseVector>
EigenLogic::diagonalizeValues(const AbstractSymmetricMatrix& symmetricMatrix) const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseSymmetricMatrix*>(&symmetricMatrix)) {
        es.compute(maybeDenseSymmetricMatrix->getDenseSymmetricMatrix(), Eigen::EigenvaluesOnly);
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseSymmetricMatrix*>(&symmetricMatrix)) {
        Eigen::MatrixXd sparseToDenseCopy = maybeSparseSymmetricMatrix->getSparseSymmetricMatrix();
        es.compute(sparseToDenseCopy, Eigen::EigenvaluesOnly);
    } else {
        throw std::bad_cast();
    }

    auto eigenvalues_ = std::make_unique<EigenDenseVector>();
    eigenvalues_->modifyDenseVector() = es.eigenvalues();

    return eigenvalues_;
}

EigenCouple
EigenLogic::diagonalizeValuesVectors(const AbstractSymmetricMatrix& symmetricMatrix) const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

    if (auto maybeDenseSymmetricMatrix =
            dynamic_cast<const EigenDenseSymmetricMatrix*>(&symmetricMatrix)) {
        es.compute(
            maybeDenseSymmetricMatrix->getDenseSymmetricMatrix(),
            Eigen::ComputeEigenvectors);
    } else if (
        auto maybeSparseSymmetricMatrix =
            dynamic_cast<const EigenSparseSymmetricMatrix*>(&symmetricMatrix)) {
        Eigen::MatrixXd sparseToDenseCopy = maybeSparseSymmetricMatrix->getSparseSymmetricMatrix();
        es.compute(sparseToDenseCopy, Eigen::ComputeEigenvectors);
    } else {
        throw std::bad_cast();
    }

    auto eigenvalues_ = std::make_unique<EigenDenseVector>();
    auto eigenvectors_ = std::make_unique<EigenDenseSemiunitaryMatrix>();
    eigenvalues_->modifyDenseVector() = es.eigenvalues();
    eigenvectors_->modifyDenseSemiunitaryMatrix() = es.eigenvectors();

    EigenCouple answer;
    answer.eigenvalues = std::move(eigenvalues_);
    answer.eigenvectors = std::move(eigenvectors_);
    return answer;
}

}  // namespace quantum::linear_algebra