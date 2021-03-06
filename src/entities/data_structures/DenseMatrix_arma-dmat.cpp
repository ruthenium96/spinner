#include <armadillo>

#include "DenseMatrix.h"

struct DenseMatrix::Impl {
    arma::dmat RawData;

  public:
    Impl() = default;
};

struct DenseVector::SubspectrumDataImpl {
    arma::dvec eigenvalues;

  public:
    SubspectrumDataImpl() = default;
};

DenseMatrix::DenseMatrix() : pImpl {std::make_unique<DenseMatrix::Impl>()} {}

DenseVector::DenseVector() : pImpl {std::make_unique<DenseVector::SubspectrumDataImpl>()} {}

DenseMatrix::~DenseMatrix() = default;
DenseMatrix::DenseMatrix(DenseMatrix&&) noexcept = default;
DenseMatrix& DenseMatrix::operator=(DenseMatrix&&) noexcept = default;
DenseVector::~DenseVector() = default;
DenseVector::DenseVector(DenseVector&&) noexcept = default;
DenseVector& DenseVector::operator=(DenseVector&&) noexcept = default;

void DenseMatrix::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->RawData(i, j) += value;
}

void DenseMatrix::assign_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->RawData(i, j) = value;
}

std::ostream& operator<<(std::ostream& os, const DenseMatrix& decomposition) {
    os << decomposition.pImpl->RawData << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const DenseVector& raw_data) {
    os << raw_data.pImpl->eigenvalues << std::endl;
    return os;
}

void DenseMatrix::resize(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize pImpl->RawData?
    pImpl->RawData.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    pImpl->RawData.fill(arma::fill::zeros);
}

void DenseMatrix::resize_with_nans(
    uint32_t matrix_in_space_basis_size_i,
    uint32_t matrix_in_space_basis_size_j) {
    // TODO: is it the fastest way to initialize pImpl->RawData?
    pImpl->RawData.resize(matrix_in_space_basis_size_i, matrix_in_space_basis_size_j);
    pImpl->RawData.fill(arma::datum::nan);
    //    for (size_t i = 0; i < pImpl->RawData.n_rows; ++i) {
    //        for (size_t j = 0; j < pImpl->RawData.n_rows; ++j) {
    //            pImpl->RawData(i, j) = NAN;
    //        }
    //    }
}

void DenseMatrix::diagonalize(DenseVector& values, DenseMatrix& vectors) const {
    arma::eig_sym(values.pImpl->eigenvalues, vectors.pImpl->RawData, pImpl->RawData);
}

void DenseMatrix::diagonalize(DenseVector& values) const {
    arma::eig_sym(values.pImpl->eigenvalues, pImpl->RawData);
}

DenseMatrix DenseMatrix::unitary_transform(const DenseMatrix& matrix_to_transform) const {
    DenseMatrix transformed_matrix;
    transformed_matrix.pImpl->RawData =
        pImpl->RawData.t() * matrix_to_transform.pImpl->RawData * pImpl->RawData;
    return std::move(transformed_matrix);
}

DenseVector DenseMatrix::return_main_diagonal() const {
    DenseVector main_diagonal;
    main_diagonal.pImpl->eigenvalues = pImpl->RawData.diag();
    return std::move(main_diagonal);
}

double DenseMatrix::operator()(uint32_t i, uint32_t j) const {
    return pImpl->RawData(i, j);
}

uint32_t DenseMatrix::size() const {
    return pImpl->RawData.n_rows;
}

uint32_t DenseMatrix::size_rows() const {
    return pImpl->RawData.n_rows;
}

uint32_t DenseMatrix::size_cols() const {
    return pImpl->RawData.n_rows;
}

uint32_t DenseVector::size() const {
    return pImpl->eigenvalues.size();
}

double DenseVector::operator()(uint32_t i) const {
    return pImpl->eigenvalues(i);
}

void DenseVector::assign_to_position(double value, uint32_t i) {
    pImpl->eigenvalues.at(i) = value;
}

void DenseVector::resize(uint32_t new_size) {
    pImpl->eigenvalues.resize(new_size);
}

std::vector<double> concatenate(const std::vector<DenseVector>& dense_vectors) {
    size_t size = 0;
    for (const auto& dense_vector : dense_vectors) {
        size += dense_vector.size();
    }
    std::vector<double> result_vector;
    result_vector.reserve(size);
    for (const auto& dense_vector : dense_vectors) {
        result_vector.insert(
            result_vector.end(),
            dense_vector.pImpl->eigenvalues.begin(),
            dense_vector.pImpl->eigenvalues.end());
    }
    return result_vector;
}

bool DenseVector::operator==(const DenseVector& rhs) const {
    return arma::approx_equal(pImpl->eigenvalues, rhs.pImpl->eigenvalues, "absdiff", 0);
}

bool DenseVector::operator!=(const DenseVector& rhs) const {
    return !(rhs == *this);
}

DenseVector DenseVector::divide_and_wise_exp(double denominator) const {
    DenseVector answer;
    answer.resize(size());
    answer.pImpl->eigenvalues = arma::exp(pImpl->eigenvalues / (denominator));
    return answer;
}

double DenseVector::dot(const DenseVector& rhs) const {
    return arma::dot(pImpl->eigenvalues, rhs.pImpl->eigenvalues);
}

DenseVector DenseVector::element_wise_multiplication(const DenseVector& rhs) const {
    DenseVector answer;
    answer.pImpl->eigenvalues = pImpl->eigenvalues % rhs.pImpl->eigenvalues;
    return answer;
}

void DenseVector::concatenate_with(const DenseVector& rhs) {
    arma::dvec tmp = arma::join_cols(pImpl->eigenvalues, rhs.pImpl->eigenvalues);
    pImpl->eigenvalues.reset();
    pImpl->eigenvalues = std::move(tmp);
}

void DenseVector::add_identical_values(size_t number, double value) {
    arma::dvec tmp;
    tmp.resize(number);
    tmp.fill(value);
    tmp = arma::join_cols(pImpl->eigenvalues, tmp);
    pImpl->eigenvalues.reset();
    pImpl->eigenvalues = std::move(tmp);
}

void DenseVector::subtract_minimum() {
    double minimum = arma::min(pImpl->eigenvalues);
    pImpl->eigenvalues -= minimum;
}