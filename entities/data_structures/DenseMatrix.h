#ifndef JULY_DENSEMATRIX_H
#define JULY_DENSEMATRIX_H

#include <iostream>
#include <memory>
#include <vector>

class DenseVector;

class DenseMatrix {
    friend DenseVector;
public:
    void add_to_position(double value, uint32_t i, uint32_t j);
    void resize(uint32_t matrix_in_space_basis_size_i, uint32_t matrix_in_space_basis_size_j);
    // TODO: is it possible to implement diagonalize() in the other place?
    void diagonalize(DenseVector& values, DenseMatrix& vectors) const;
    void diagonalize(DenseVector& values) const;
    [[nodiscard]] DenseVector return_main_diagonal() const;
    [[nodiscard]] DenseMatrix unitary_transform(const DenseMatrix& matrix_to_transform) const;

    [[nodiscard]] uint32_t size() const;
    double operator()(uint32_t i, uint32_t j) const;

    DenseMatrix();
    DenseMatrix(const DenseMatrix&) = delete;
    DenseMatrix& operator=(const DenseMatrix&) = delete;
    DenseMatrix(DenseMatrix&&) noexcept;
    DenseMatrix& operator=(DenseMatrix&&) noexcept;
    ~DenseMatrix();

    friend std::ostream &operator<<(std::ostream &os, const DenseMatrix &raw_data);

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};

class DenseVector {
    friend DenseMatrix;
public:
    DenseVector();
    DenseVector(const DenseVector&) = delete;
    DenseVector& operator=(const DenseVector&) = delete;
    DenseVector(DenseVector&&) noexcept;
    DenseVector& operator=(DenseVector&&) noexcept;
    ~DenseVector();

    [[nodiscard]] uint32_t size() const;

    friend std::ostream &operator<<(std::ostream &os, const DenseVector &raw_data);

    // TODO: it is a temporary solution, fix it
    friend std::vector<double> concatenate(const std::vector<DenseVector>&);
private:
    class SubspectrumDataImpl;
    std::unique_ptr<SubspectrumDataImpl> pImpl;
};

#endif //JULY_DENSEMATRIX_H
