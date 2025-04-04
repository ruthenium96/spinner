#include "ArmaDenseVector.h"
#include <stdexcept>

namespace quantum::linear_algebra {

template <typename T>
void ArmaDenseVector<T>::resize(uint32_t new_size) {
    vector_.resize(new_size);
}

template <typename T>
void ArmaDenseVector<T>::makeRandomUnitVector(uint32_t size) {
    if (!vector_.empty()) {
        throw std::invalid_argument("Trying to make random unit vector from non-empty vector.");
    }
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution d{0.0, 1.0};

    vector_.resize(size);
    for (uint32_t i = 0; i < size; ++i) {
        vector_(i) = d(gen);
    }

    double vec_norm = arma::vecnorm(vector_);
    vector_ /= vec_norm;
}

template <typename T>
void ArmaDenseVector<T>::concatenate_with(const std::unique_ptr<AbstractDenseVector>& rhs) {
    auto rhs_ = downcast_ptr(rhs);

    arma::Col<T> tmp = arma::join_cols(vector_, rhs_->vector_);
    vector_.reset();
    vector_ = std::move(tmp);
}

template <typename T>
void ArmaDenseVector<T>::add_identical_values(size_t number, double value) {
    arma::Col<T> tmp;
    tmp.resize(number);
    tmp.fill(value);
    tmp = arma::join_cols(vector_, tmp);
    vector_.reset();
    vector_ = std::move(tmp);
}

template <typename T>
void ArmaDenseVector<T>::subtract_minimum() {
    double minimum = arma::min(vector_);
    vector_ -= minimum;
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaDenseVector<T>::multiply_by(double multiplier) const {
    auto answer = std::make_unique<ArmaDenseVector<T>>();
    answer->resize(size());
    answer->vector_ = vector_ * multiplier;
    return answer;
}

template <typename T>
void ArmaDenseVector<T>::wise_exp() {
    vector_ = arma::exp(vector_);
}

template <typename T>
double ArmaDenseVector<T>::dot(const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    return arma::dot(vector_, rhs_->vector_);
}

template <typename T>
double ArmaDenseVector<T>::triple_dot(
    const std::unique_ptr<AbstractDenseVector>& second, 
    const std::unique_ptr<AbstractDenseVector>& third) const {
    auto second_ = downcast_ptr(second);
    auto third_ = downcast_ptr(third);

    return arma::dot(vector_, second_->vector_ % third_->vector_);
}

template <typename T>
std::unique_ptr<AbstractDenseVector> ArmaDenseVector<T>::element_wise_multiplication(
    const std::unique_ptr<AbstractDenseVector>& rhs) const {
    auto rhs_ = downcast_ptr(rhs);

    auto answer = std::make_unique<ArmaDenseVector>();
    answer->vector_ = vector_ % rhs_->vector_;
    return answer;
}

template <typename T>
uint32_t ArmaDenseVector<T>::size() const {
    return vector_.size();
}

template <typename T>
double ArmaDenseVector<T>::at(uint32_t i) const {
    return vector_.at(i);
}

template <typename T>
const ArmaDenseVector<T>*
ArmaDenseVector<T>::downcast_ptr(const std::unique_ptr<AbstractDenseVector>& ptr) {
    auto answer = dynamic_cast<const ArmaDenseVector<T>*>(ptr.get());
    if (answer == nullptr) {
        throw std::bad_cast();
    }
    return answer;
}

template <typename T>
void ArmaDenseVector<T>::print(std::ostream& os) const {
    os << vector_ << std::endl;
}

template <typename T>
arma::Col<T>& ArmaDenseVector<T>::modifyDenseVector() {
    return vector_;
}

template <typename T>
const arma::Col<T>& ArmaDenseVector<T>::getDenseVector() const {
    return vector_;
}

template class ArmaDenseVector<double>;
template class ArmaDenseVector<float>;
}  // namespace quantum::linear_algebra