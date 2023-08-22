#ifndef SPINNER_MULTIPLICITYDIRECTSUM_H
#define SPINNER_MULTIPLICITYDIRECTSUM_H

#include <vector>

#include "Multiplicity.h"

namespace spin_algebra {

class MultiplicityDirectSum {
  private:
    std::vector<Multiplicity> multiplicities_;

  public:
    MultiplicityDirectSum() = default;
    explicit MultiplicityDirectSum(Multiplicity single_multiplicity);
    explicit MultiplicityDirectSum(std::vector<Multiplicity> multiplicities);
    MultiplicityDirectSum(std::initializer_list<Multiplicity> l);

    const std::vector<Multiplicity>& getMultiplicities() const;

    MultiplicityDirectSum operator+(const MultiplicityDirectSum& bs) const;
    MultiplicityDirectSum& operator+=(const MultiplicityDirectSum& bs);
    MultiplicityDirectSum operator*(const MultiplicityDirectSum& bs) const;
    MultiplicityDirectSum& operator*=(const MultiplicityDirectSum& bs);

    bool operator==(const MultiplicityDirectSum& rhs) const;
    bool operator!=(const MultiplicityDirectSum& rhs) const;
};

}  // namespace spin_algebra
#endif  //SPINNER_MULTIPLICITYDIRECTSUM_H
