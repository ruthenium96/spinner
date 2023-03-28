#include "MultiplicityDirectSum.h"

#include <algorithm>
#include <utility>

namespace spin_algebra {

MultiplicityDirectSum::MultiplicityDirectSum(Multiplicity single_multiplicity) {
    multiplicities_.push_back(single_multiplicity);
}

MultiplicityDirectSum::MultiplicityDirectSum(std::vector<Multiplicity> multiplicities) {
    multiplicities_ = std::move(multiplicities);
}

MultiplicityDirectSum::MultiplicityDirectSum(std::initializer_list<Multiplicity> l) :
    multiplicities_(l) {}

const std::vector<Multiplicity>& MultiplicityDirectSum::getMultiplicities() const {
    return multiplicities_;
}

MultiplicityDirectSum MultiplicityDirectSum::operator+(const MultiplicityDirectSum& bs) const {
    MultiplicityDirectSum answer;
    answer.multiplicities_.reserve(this->multiplicities_.size() + bs.multiplicities_.size());
    copy(
        this->multiplicities_.begin(),
        this->multiplicities_.end(),
        back_inserter(answer.multiplicities_));
    answer += bs;
    return answer;
}

MultiplicityDirectSum MultiplicityDirectSum::operator*(const MultiplicityDirectSum& bs) const {
    MultiplicityDirectSum answer;

    for (const auto& a : this->multiplicities_) {
        for (const auto& b : bs.multiplicities_) {
            Multiplicity min_multiplicity = std::min(a, b);
            Multiplicity max_multiplicity = std::max(a, b);
            for (Multiplicity new_multiplicity = max_multiplicity - min_multiplicity + 1;
                 new_multiplicity < max_multiplicity + min_multiplicity;
                 new_multiplicity += 2) {
                answer.multiplicities_.push_back(new_multiplicity);
            }
        }
    }

    return answer;
}

MultiplicityDirectSum& MultiplicityDirectSum::operator+=(const MultiplicityDirectSum& bs) {
    this->multiplicities_.reserve(this->multiplicities_.size() + bs.multiplicities_.size());
    copy(
        bs.multiplicities_.begin(),
        bs.multiplicities_.end(),
        back_inserter(this->multiplicities_));
    return *this;
}

MultiplicityDirectSum& MultiplicityDirectSum::operator*=(const MultiplicityDirectSum& bs) {
    std::vector<Multiplicity> answer_vector;

    for (const auto& a : this->multiplicities_) {
        for (const auto& b : bs.multiplicities_) {
            Multiplicity min_multiplicity = std::min(a, b);
            Multiplicity max_multiplicity = std::max(a, b);
            for (uint16_t new_multiplicity = max_multiplicity - min_multiplicity + 1;
                 new_multiplicity < max_multiplicity + min_multiplicity;
                 new_multiplicity += 2) {
                answer_vector.push_back(new_multiplicity);
            }
        }
    }

    this->multiplicities_ = std::move(answer_vector);

    return *this;
}

bool MultiplicityDirectSum::operator==(const MultiplicityDirectSum& rhs) const {
    std::vector<Multiplicity> this_sorted = this->multiplicities_;
    std::vector<Multiplicity> rhs_sorted = rhs.multiplicities_;

    std::sort(this_sorted.begin(), this_sorted.end());
    std::sort(rhs_sorted.begin(), rhs_sorted.end());

    return this_sorted == rhs_sorted;
}

bool MultiplicityDirectSum::operator!=(const MultiplicityDirectSum& rhs) const {
    return !(rhs == *this);
}

}  // namespace spin_algebra