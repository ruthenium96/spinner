#ifndef SPINNER_TERM_H
#define SPINNER_TERM_H

#include <cstdint>
#include <memory>

#include "src/common/lexicographic/IndexConverter.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"
#include "src/model/NumericalParameters.h"

namespace model::operators {
class Term {
  public:
    // This method is required for deep copy of std::vector<std::unique_ptr<Term>>
    virtual std::unique_ptr<Term> clone() const = 0;
    virtual void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        uint32_t index_of_vector) const = 0;
    virtual ~Term() = default;
};

class ZeroCenterTerm : public Term {
  public:
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        uint32_t index_of_vector) const override = 0;
    ~ZeroCenterTerm() override = default;
};

class OneCenterTerm : public Term {
  private:
    size_t numberOfCenters_;
  public:
    size_t getNumberOfCenters() const {
        return numberOfCenters_;
    }
    explicit OneCenterTerm(size_t numberOfCenters) : numberOfCenters_(numberOfCenters) {}
    virtual void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a) const = 0;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector) const override {
        for (int center_a = 0; center_a < getNumberOfCenters(); ++center_a) {
            construct(
                matrix_in_lexicografical_basis,
                index_of_vector,
                center_a);
        }
    };
    ~OneCenterTerm() override = default;
};

class TwoCenterTerm : public Term {
  private:
    size_t numberOfCenters_;
  public:
    size_t getNumberOfCenters() const {
        return numberOfCenters_;
    }
    explicit TwoCenterTerm(size_t numberOfCenters) : numberOfCenters_(numberOfCenters) {}
    virtual void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix&
            matrix_in_lexicografical_basis,
        uint32_t index_of_vector,
        uint32_t center_a,
        uint32_t center_b) const = 0;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix& matrix_in_lexicografical_basis,
        uint32_t index_of_vector) const override {
        for (int center_a = 0; center_a < getNumberOfCenters(); ++center_a) {
            for (int center_b = center_a + 1; center_b < getNumberOfCenters(); ++center_b) {
                construct(
                    matrix_in_lexicografical_basis,
                    index_of_vector,
                    center_a,
                    center_b);
            }
        }
    };
    ~TwoCenterTerm() override = default;
};
}  // namespace model::operators

#endif  //SPINNER_TERM_H
