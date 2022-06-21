#ifndef SPINNER_ARMAVECTOR_H
#define SPINNER_ARMAVECTOR_H

#include <armadillo>
#include <cstddef>
#include <cstdint>

#include "src/entities/data_structures/AbstractVector.h"

namespace quantum::linear_algebra {
class ArmaFactory;
class ArmaMatrix;
class ArmaVector: public AbstractVector {
    friend ArmaFactory;
    friend ArmaMatrix;

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
    bool operator==(const std::unique_ptr<AbstractVector>& rhs) const override;
    bool operator!=(const std::unique_ptr<AbstractVector>& rhs) const override;
    void print(std::ostream& os) const override;

  private:
    // c-like pointers are necessary to avoid double-free error
    static const ArmaVector* downcast_ptr(const std::unique_ptr<AbstractVector>& ptr);
    static const ArmaVector* downcast_ptr(std::unique_ptr<const AbstractVector>& ptr);
    static ArmaVector* downcast_ptr(std::unique_ptr<AbstractVector>& ptr);
    arma::dvec vector_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_ARMAVECTOR_H
