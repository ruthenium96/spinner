#include "EigenVector.h"

namespace quantum::linear_algebra {

void EigenVector::assign_to_position(double value, uint32_t i) {
    vector_(i) = value;
}

void EigenVector::resize(uint32_t new_size) {
    vector_.resize(new_size);
}

void EigenVector::concatenate_with(const std::unique_ptr<AbstractVector>& rhs) {
    auto rhs_ = downcast_ptr(rhs);

    Eigen::VectorXd tmp(vector_.size() + rhs_->vector_.size());
    tmp << vector_, rhs_->vector_;

    vector_ = tmp;
}

void EigenVector::add_identical_values(size_t number, double value) {
    Eigen::VectorXd tmp(number);
    tmp.fill(value);

    Eigen::VectorXd tmptmp(vector_.size() + number);
    tmptmp << vector_, tmp;
    vector_ = tmptmp;
}

void EigenVector::subtract_minimum() {
    double minimum = vector_.minCoeff();
    vector_.array() -= minimum;
}

std::unique_ptr<AbstractVector> EigenVector::divide_and_wise_exp(double denominator) const {
    auto answer = std::make_unique<EigenVector>();
    answer->resize(vector_.size());

    answer->vector_.array() = exp(vector_.array() / denominator);

    return answer;
}

double EigenVector::dot(const std::unique_ptr<AbstractVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return vector_.dot(rhs_->vector_);
}

std::unique_ptr<AbstractVector>
EigenVector::element_wise_multiplication(const std::unique_ptr<AbstractVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    auto answer = std::make_unique<EigenVector>();
    answer->resize(vector_.size());

    answer->vector_.array() = vector_.array() * rhs_->vector_.array();

    return answer;
}

uint32_t EigenVector::size() const {
    return vector_.size();
}

double EigenVector::at(uint32_t i) const {
    return vector_(i);
}

void EigenVector::print(std::ostream& os) const {
    os << vector_ << std::endl;
}

bool EigenVector::operator==(const std::unique_ptr<AbstractVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return vector_ == rhs_->vector_;
}

bool EigenVector::operator!=(const std::unique_ptr<AbstractVector>& rhs) const {
    return !(*this == rhs);
}

const EigenVector* EigenVector::downcast_ptr(const std::unique_ptr<AbstractVector>& ptr) {
    auto answer = dynamic_cast<const EigenVector*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

}  // namespace quantum::linear_algebra