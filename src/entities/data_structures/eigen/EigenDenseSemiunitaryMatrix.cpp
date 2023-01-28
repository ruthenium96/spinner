#include "EigenDenseSemiunitaryMatrix.h"

#include "EigenDenseSymmetricMatrix.h"
#include "EigenDenseVector.h"

namespace quantum::linear_algebra {

uint32_t EigenDenseSemiunitaryMatrix::size_rows() const {
    return denseSemiunitaryMatrix_.rows();
}

uint32_t EigenDenseSemiunitaryMatrix::size_cols() const {
    return denseSemiunitaryMatrix_.cols();
}

void EigenDenseSemiunitaryMatrix::resize(size_t size_rows, size_t size_cols) {
    denseSemiunitaryMatrix_.resize(size_rows, size_cols);
    denseSemiunitaryMatrix_.fill(0);
}

double EigenDenseSemiunitaryMatrix::at(uint32_t i, uint32_t j) const {
    return denseSemiunitaryMatrix_(i, j);
}

void EigenDenseSemiunitaryMatrix::print(std::ostream& os) const {
    os << denseSemiunitaryMatrix_ << std::endl;
}

std::unique_ptr<AbstractDenseVector>
EigenDenseSemiunitaryMatrix::unitaryTransformAndReturnMainDiagonal(
    const std::unique_ptr<AbstractSymmetricMatrix>& matrix_to_transform) const {
    auto answer = std::make_unique<EigenDenseVector>();

    auto matrix_to_transform_dense =
        dynamic_cast<const EigenDenseSymmetricMatrix*>(matrix_to_transform.get());
    if (matrix_to_transform_dense != nullptr) {
        Eigen::MatrixXd transformed_matrix = denseSemiunitaryMatrix_.transpose()
            * matrix_to_transform_dense->getDenseSymmetricMatrix().transpose()
            * denseSemiunitaryMatrix_;

        answer->modifyDenseVector() = transformed_matrix.diagonal();
        return answer;
    }
    throw std::bad_cast();
}

const Eigen::MatrixXd& EigenDenseSemiunitaryMatrix::getDenseSemiunitaryMatrix() const {
    return denseSemiunitaryMatrix_;
}

Eigen::MatrixXd& EigenDenseSemiunitaryMatrix::modifyDenseSemiunitaryMatrix() {
    return denseSemiunitaryMatrix_;
}
}  // namespace quantum::linear_algebra