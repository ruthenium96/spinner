#include "ArmaDenseVector.h"

namespace quantum::linear_algebra {

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

std::unique_ptr<AbstractDenseVector> ArmaDenseVector::multiply_by(double multiplier) const {
    auto answer = std::make_unique<ArmaDenseVector>();
    answer->resize(size());
    answer->vector_ = vector_ * multiplier;
    return answer;
}

void ArmaDenseVector::wise_exp() {
    vector_ = arma::exp(vector_);
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

arma::dvec& ArmaDenseVector::modifyDenseVector() {
    return vector_;
}

const arma::dvec& ArmaDenseVector::getDenseVector() {
    return vector_;
}

}  // namespace quantum::linear_algebra