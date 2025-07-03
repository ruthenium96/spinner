#ifndef SPINNER_UNCERTAINVALUE_H
#define SPINNER_UNCERTAINVALUE_H

#include <array>
#include <cmath>

namespace common {

enum UncertaintySources {
    FTLM,
    EXPERIMENT,
    FIT // TODO: it is not real source of uncertainty, it is just mix of previous two.
        // we need to both modify loss function and use Delta Method. 
};

class UncertainValue {
  private:
    double mean_;
    static constexpr size_t NUMBER_OF_SOURCES = 3;
    std::array<double, NUMBER_OF_SOURCES> stdevs_;
    // while FTLM uncertainties are correlated, EXPERIMENT -- not:
    std::array<double, NUMBER_OF_SOURCES> correlation_coeffs_ = {1.0, 0.0, 0.0};
  public:
    UncertainValue();
    explicit UncertainValue(double mean);
    UncertainValue(double mean, double stdev, UncertaintySources source);
    UncertainValue(double mean, std::array<double, NUMBER_OF_SOURCES> stdevs);

    double mean() const;
    const std::array<double, NUMBER_OF_SOURCES>& stdevs() const;
    double stdev_total() const;

    UncertainValue operator-() const;

    UncertainValue operator+(const UncertainValue& rhs) const;
    UncertainValue& operator+=(const UncertainValue& rhs);

    UncertainValue operator-(const UncertainValue& rhs) const;
    UncertainValue& operator-=(const UncertainValue& rhs);

    UncertainValue operator*(const UncertainValue& rhs) const;
    UncertainValue& operator*=(const UncertainValue& rhs);
    UncertainValue operator*(double rhs) const;
    UncertainValue& operator*=(double rhs);
    friend UncertainValue operator*(double lhs, const UncertainValue& rhs);

    UncertainValue operator/(const UncertainValue& rhs) const;
    UncertainValue& operator/=(const UncertainValue& rhs);
    UncertainValue operator/(double rhs) const;
    UncertainValue& operator/=(double rhs);

    static UncertainValue sqrt(const UncertainValue& value);
    static UncertainValue inv(const UncertainValue& value);
};

} // common

#endif //SPINNER_UNCERTAINVALUE_H
