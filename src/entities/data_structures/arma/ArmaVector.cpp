#include "ArmaVector.h"

namespace quantum::linear_algebra {

void ArmaVector::assign_to_position(double value, uint32_t i) {
    vector_.at(i) = value;
}

void ArmaVector::resize(uint32_t new_size) {
    vector_.resize(new_size);
}

void ArmaVector::concatenate_with(const std::unique_ptr<AbstractVector>& rhs) {
    auto rhs_ = downcast_ptr(rhs);

    arma::dvec tmp = arma::join_cols(vector_, rhs_->vector_);
    vector_.reset();
    vector_ = std::move(tmp);
}

void ArmaVector::add_identical_values(size_t number, double value) {
    arma::dvec tmp;
    tmp.resize(number);
    tmp.fill(value);
    tmp = arma::join_cols(vector_, tmp);
    vector_.reset();
    vector_ = std::move(tmp);
}

void ArmaVector::subtract_minimum() {
    double minimum = arma::min(vector_);
    vector_ -= minimum;
}

std::unique_ptr<AbstractVector> ArmaVector::divide_and_wise_exp(double denominator) const {
    auto answer = std::make_unique<ArmaVector>();
    answer->resize(size());
    answer->vector_ = arma::exp(vector_ / (denominator));
    return answer;
}

double ArmaVector::dot(const std::unique_ptr<AbstractVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return arma::dot(vector_, rhs_->vector_);
}

std::unique_ptr<AbstractVector>
ArmaVector::element_wise_multiplication(const std::unique_ptr<AbstractVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    auto answer = std::make_unique<ArmaVector>();
    answer->vector_ = vector_ % rhs_->vector_;
    return answer;
}

uint32_t ArmaVector::size() const {
    return vector_.size();
}

double ArmaVector::at(uint32_t i) const {
    return vector_.at(i);
}

bool ArmaVector::operator==(const std::unique_ptr<AbstractVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return arma::approx_equal(vector_, rhs_->vector_, "absdiff", 0);
}

bool ArmaVector::operator!=(const std::unique_ptr<AbstractVector>& rhs) const {
    return !(*this == rhs);
}

const ArmaVector* ArmaVector::downcast_ptr(const std::unique_ptr<AbstractVector>& ptr) {
    return dynamic_cast<const ArmaVector*>(ptr.get());
}

void ArmaVector::print(std::ostream& os) const {
    os << vector_ << std::endl;
}

}  // namespace quantum::linear_algebra