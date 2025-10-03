#ifndef SPINNER_ABSTRACTDENSEVECTOR_H
#define SPINNER_ABSTRACTDENSEVECTOR_H

#include <cstdint>
#include <memory>

namespace quantum::linear_algebra {
class AbstractDenseVector {
  public:
    virtual void concatenate_with(const std::unique_ptr<AbstractDenseVector>& rhs) = 0;
    virtual void add_identical_values(size_t number, double value) = 0;
    virtual void subtract_minimum() = 0;
    virtual void wise_exp() = 0;

    virtual std::unique_ptr<AbstractDenseVector> multiply_by(double multiplier) const = 0;
    virtual double dot(const std::unique_ptr<AbstractDenseVector>& rhs) const = 0;
    // einsum of "iii->" type:
    virtual double triple_dot(const std::unique_ptr<AbstractDenseVector>& second, 
      const std::unique_ptr<AbstractDenseVector>& third) const = 0;
    virtual std::unique_ptr<AbstractDenseVector>
    element_wise_multiplication(const std::unique_ptr<AbstractDenseVector>& rhs) const = 0;

    virtual uint32_t size() const = 0;
    virtual double at(uint32_t i) const = 0;

    virtual void print(std::ostream& os) const = 0;

    virtual ~AbstractDenseVector() = default;
};
}  // namespace quantum::linear_algebra

#endif  //SPINNER_ABSTRACTDENSEVECTOR_H