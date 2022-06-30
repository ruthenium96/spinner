#ifndef SPINNER_EIGENVECTOR_H
#define SPINNER_EIGENVECTOR_H

#include <Eigen/Dense>

#include "src/entities/data_structures/AbstractVector.h"

namespace quantum::linear_algebra {
class EigenFactory;
class EigenMatrix;
class EigenVector: public AbstractVector {
    friend EigenFactory;
    friend EigenMatrix;

  public:
    void assign_to_position(double value, uint32_t i) override;
    void resize(uint32_t new_size) override;
    void concatenate_with(const std::unique_ptr<AbstractVector>& rhs) override;
    void add_identical_values(size_t number, double value) override;
    void subtract_minimum() override;
    std::unique_ptr<AbstractVector> divide_and_wise_exp(double denominator) const override;
    double dot(const std::unique_ptr<AbstractVector>& rhs) const override;
    std::unique_ptr<AbstractVector>
    element_wise_multiplication(const std::unique_ptr<AbstractVector>& rhs) const override;
    uint32_t size() const override;
    double at(uint32_t i) const override;
    void print(std::ostream& os) const override;

    bool operator==(const std::unique_ptr<AbstractVector>& rhs) const override;
    bool operator!=(const std::unique_ptr<AbstractVector>& rhs) const override;

  private:
    Eigen::VectorXd vector_;

    // c-like pointers are necessary to avoid double-free error
    static const EigenVector* downcast_ptr(const std::unique_ptr<AbstractVector>& ptr);
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EIGENVECTOR_H
