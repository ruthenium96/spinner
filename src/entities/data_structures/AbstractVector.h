#ifndef SPINNER_ABSTRACTVECTOR_H
#define SPINNER_ABSTRACTVECTOR_H

#include <cstdint>
#include <memory>
#include <vector>

namespace quantum::linear_algebra {
class AbstractMatrix;
class AbstractVector {
  public:
    virtual void assign_to_position(double value, uint32_t i) = 0;
    virtual void resize(uint32_t new_size) = 0;

    virtual void concatenate_with(const std::unique_ptr<AbstractVector>& rhs) = 0;
    virtual void add_identical_values(size_t number, double value) = 0;
    virtual void subtract_minimum() = 0;

    virtual std::unique_ptr<AbstractVector> divide_and_wise_exp(double denominator) const = 0;
    virtual double dot(const std::unique_ptr<AbstractVector>& rhs) const = 0;
    virtual std::unique_ptr<AbstractVector>
    element_wise_multiplication(const std::unique_ptr<AbstractVector>& rhs) const = 0;

    virtual uint32_t size() const = 0;
    virtual double at(uint32_t i) const = 0;

    virtual void print(std::ostream& os) const = 0;

    virtual bool operator==(const std::unique_ptr<AbstractVector>& rhs) const = 0;
    virtual bool operator!=(const std::unique_ptr<AbstractVector>& rhs) const = 0;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTVECTOR_H