#include "EigenDenseVector.h"

namespace quantum::linear_algebra {

void EigenDenseVector::resize(uint32_t new_size) {
    vector_.resize(new_size);
}

void EigenDenseVector::concatenate_with(const std::unique_ptr<AbstractDenseVector>& rhs) {
    auto rhs_ = downcast_ptr(rhs);

    Eigen::VectorXd tmp(vector_.size() + rhs_->vector_.size());
    tmp << vector_, rhs_->vector_;

    vector_ = tmp;
}

void EigenDenseVector::add_identical_values(size_t number, double value) {
    Eigen::VectorXd tmp(number);
    tmp.fill(value);

    Eigen::VectorXd tmptmp(vector_.size() + number);
    tmptmp << vector_, tmp;
    vector_ = tmptmp;
}

void EigenDenseVector::subtract_minimum() {
    double minimum = vector_.minCoeff();
    vector_.array() -= minimum;
}

std::unique_ptr<AbstractDenseVector> EigenDenseVector::multiply_by(double multiplier) const {
    auto answer = std::make_unique<EigenDenseVector>();
    answer->resize(vector_.size());

    answer->vector_.array() = vector_.array() * multiplier;

    return answer;
}

void EigenDenseVector::wise_exp() {
    vector_.array() = exp(vector_.array());
}

double EigenDenseVector::dot(const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return vector_.dot(rhs_->vector_);
}

std::unique_ptr<AbstractDenseVector> EigenDenseVector::element_wise_multiplication(
    const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    auto answer = std::make_unique<EigenDenseVector>();
    answer->resize(vector_.size());

    answer->vector_.array() = vector_.array() * rhs_->vector_.array();

    return answer;
}

uint32_t EigenDenseVector::size() const {
    return vector_.size();
}

double EigenDenseVector::at(uint32_t i) const {
    return vector_(i);
}

void EigenDenseVector::print(std::ostream& os) const {
    os << vector_ << std::endl;
}

Eigen::VectorXd& EigenDenseVector::modifyDenseVector() {
    return vector_;
}

const Eigen::VectorXd& EigenDenseVector::getDenseVector() {
    return vector_;
}

const EigenDenseVector*
EigenDenseVector::downcast_ptr(const std::unique_ptr<AbstractDenseVector>& ptr) {
    auto answer = dynamic_cast<const EigenDenseVector*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

}  // namespace quantum::linear_algebra