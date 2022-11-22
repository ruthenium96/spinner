#include "ArmaDenseVector.h"

namespace quantum::linear_algebra {

void ArmaDenseVector::assign_to_position(double value, uint32_t i) {
    vector_.at(i) = value;
}

void ArmaDenseVector::resize(uint32_t new_size) {
    vector_.resize(new_size);
}

void ArmaDenseVector::concatenate_with(const std::unique_ptr<AbstractDenseVector>& rhs) {
    auto rhs_ = downcast_ptr(rhs);

    arma::dvec tmp = arma::join_cols(vector_, rhs_->vector_);
    vector_.reset();
    vector_ = std::move(tmp);
}

void ArmaDenseVector::add_identical_values(size_t number, double value) {
    arma::dvec tmp;
    tmp.resize(number);
    tmp.fill(value);
    tmp = arma::join_cols(vector_, tmp);
    vector_.reset();
    vector_ = std::move(tmp);
}

void ArmaDenseVector::subtract_minimum() {
    double minimum = arma::min(vector_);
    vector_ -= minimum;
}

std::unique_ptr<AbstractDenseVector>
ArmaDenseVector::divide_and_wise_exp(double denominator) const {
    auto answer = std::make_unique<ArmaDenseVector>();
    answer->resize(size());
    answer->vector_ = arma::exp(vector_ / (denominator));
    return answer;
}

double ArmaDenseVector::dot(const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return arma::dot(vector_, rhs_->vector_);
}

std::unique_ptr<AbstractDenseVector> ArmaDenseVector::element_wise_multiplication(
    const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    auto answer = std::make_unique<ArmaDenseVector>();
    answer->vector_ = vector_ % rhs_->vector_;
    return answer;
}

uint32_t ArmaDenseVector::size() const {
    return vector_.size();
}

double ArmaDenseVector::at(uint32_t i) const {
    return vector_.at(i);
}

bool ArmaDenseVector::operator==(const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return arma::approx_equal(vector_, rhs_->vector_, "absdiff", 0);
}

bool ArmaDenseVector::operator!=(const std::unique_ptr<AbstractDenseVector>& rhs) const {
    return !(*this == rhs);
}

const ArmaDenseVector*
ArmaDenseVector::downcast_ptr(const std::unique_ptr<AbstractDenseVector>& ptr) {
    auto answer = dynamic_cast<const ArmaDenseVector*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

void ArmaDenseVector::print(std::ostream& os) const {
    os << vector_ << std::endl;
}

}  // namespace quantum::linear_algebra