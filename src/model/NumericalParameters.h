#ifndef SPINNER_NUMERICALPARAMETERS_H
#define SPINNER_NUMERICALPARAMETERS_H

#include <vector>

namespace model {
template<class T>
class TwoDNumericalParameters {
  private:
    std::vector<std::vector<T>> data_;

  public:
    explicit TwoDNumericalParameters(size_t size, T default_value) {
        data_ = std::vector<std::vector<T>>(size, std::vector<T>(size, default_value));
    }
    size_t size() const {
        return data_.size();
    }
    T& at(size_t i, size_t j) {
        return data_.at(i).at(j);
    }
    const T& at(size_t i, size_t j) const {
        return data_.at(i).at(j);
    }
};

template<class T>
class OneDNumericalParameters {
  private:
    std::vector<T> data_;

  public:
    explicit OneDNumericalParameters(size_t size, T default_value) {
        data_ = std::vector<T>(size, default_value);
    }
    size_t size() const {
        return data_.size();
    }
    T& at(size_t i) {
        return data_.at(i);
    }
    const T& at(size_t i) const {
        return data_.at(i);
    }
};
}  // namespace model

#endif  //SPINNER_NUMERICALPARAMETERS_H
