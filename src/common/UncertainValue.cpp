#include "UncertainValue.h"
#include <stdexcept>

namespace {
template <typename T> 
inline constexpr int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
}

namespace common {

UncertainValue::UncertainValue() : mean_(0.0), stdevs_({}) {}

UncertainValue::UncertainValue(double mean) : mean_(mean), stdevs_({}) {}

UncertainValue::UncertainValue(double mean, double stdev, UncertaintySources source) : mean_(mean), stdevs_({}) {
    if (stdev < 0) {
        throw std::invalid_argument("stdev cannot be negative!");
    }
    stdevs_[source] = stdev;
}

UncertainValue::UncertainValue(double mean, std::array<double, NUMBER_OF_SOURCES> stdevs) {
    mean_ = mean;
    for (double stdev : stdevs) {
        if (stdev < 0) {
            throw std::invalid_argument("stdev cannot be negative!");
        }
    }
    stdevs_ = std::move(stdevs);
}

double UncertainValue::mean() const {
    return mean_;
}

const std::array<double, UncertainValue::NUMBER_OF_SOURCES>& UncertainValue::stdevs() const {
    return stdevs_;
}

double UncertainValue::stdev_total() const {
    double sd_squared = 0.0;
    for (const auto el : stdevs_) {
        sd_squared += el * el;
    }
    return std::sqrt(sd_squared);
}

UncertainValue UncertainValue::operator-() const {
    UncertainValue answer = *this;
    answer *= -1.0;
    return answer;
}

UncertainValue UncertainValue::operator+(const UncertainValue& rhs) const {
    UncertainValue answer = *this;
    answer += rhs;
    return answer;
}

UncertainValue& UncertainValue::operator+=(const UncertainValue& rhs) {
    double new_mean = mean_ + rhs.mean_;
    std::array<double, NUMBER_OF_SOURCES> new_stdevs = {};
    for (size_t i = 0; i < NUMBER_OF_SOURCES; ++i) {
        new_stdevs[i] = std::sqrt(
            stdevs_[i] * stdevs_[i] + 
            rhs.stdevs_[i] * rhs.stdevs_[i] +
            2.0 * stdevs_[i] * rhs.stdevs_[i] * correlation_coeffs_[i] * rhs.correlation_coeffs_[i]);

        correlation_coeffs_[i] = sgn(correlation_coeffs_[i] * stdevs_[i] + 
            rhs.correlation_coeffs_[i] * rhs.stdevs_[i]);
    }
    mean_ = new_mean;
    stdevs_ = new_stdevs;
    return *this;
}

UncertainValue UncertainValue::operator-(const UncertainValue& rhs) const {
    UncertainValue answer = *this;
    answer -= rhs;
    return answer;
}

UncertainValue& UncertainValue::operator-=(const UncertainValue& rhs) {
    *this += -rhs;
    return *this;
}

UncertainValue UncertainValue::operator*(const UncertainValue& rhs) const {
    UncertainValue answer = *this;
    answer *= rhs;
    return answer;
}

UncertainValue& UncertainValue::operator*=(const UncertainValue& rhs) {
    double new_mean = mean_ * rhs.mean_;
    std::array<double, NUMBER_OF_SOURCES> new_stdevs = {};
    for (size_t i = 0; i < NUMBER_OF_SOURCES; ++i) {
        double first_factor = stdevs_[i] * std::abs(rhs.mean_);
        double second_factor = rhs.stdevs_[i] * std::abs(mean_);
        new_stdevs[i] = std::sqrt(
            first_factor * first_factor +
            second_factor * second_factor +
            2 * first_factor * second_factor * correlation_coeffs_[i] * rhs.correlation_coeffs_[i]
        );
        correlation_coeffs_[i] = sgn(correlation_coeffs_[i] * stdevs_[i] * rhs.mean_ + 
            rhs.correlation_coeffs_[i] * rhs.stdevs_[i] * mean_);
    }
    mean_ = new_mean;
    stdevs_ = new_stdevs;
    return *this;
}

UncertainValue operator*(double lhs, const UncertainValue& rhs) {
    UncertainValue answer = rhs;
    answer *= lhs;
    return answer;
}

UncertainValue UncertainValue::operator/(const UncertainValue& rhs) const {
    UncertainValue answer = *this;
    answer /= rhs;
    return answer;
}

UncertainValue& UncertainValue::operator/=(const UncertainValue& rhs) {
    *this *= inv(rhs);
    return *this;
}

UncertainValue UncertainValue::sqrt(const UncertainValue& value) {
    UncertainValue answer;
    answer.mean_ = std::sqrt(value.mean_);
    for (size_t i = 0; i < NUMBER_OF_SOURCES; ++i) {
        answer.stdevs_[i] = value.stdevs_[i] / (2 * std::sqrt(value.mean_));
    }
    answer.correlation_coeffs_ = value.correlation_coeffs_;
    return answer;
}

UncertainValue UncertainValue::inv(const UncertainValue& value) {
    UncertainValue answer;
    answer.mean_ = 1 / value.mean_;
    for (size_t i = 0; i < NUMBER_OF_SOURCES; ++i) {
        answer.stdevs_[i] = value.stdevs_[i] / (value.mean_ * value.mean_);
        answer.correlation_coeffs_[i] = -1.0 * value.correlation_coeffs_[i];
    }
    return answer;
}

} // common