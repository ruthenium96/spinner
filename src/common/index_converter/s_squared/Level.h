#ifndef SPINNER_S2LEVEL_H
#define SPINNER_S2LEVEL_H

#include "src/spin_algebra/Multiplicity.h"

#include <compare>
#include <memory>
#include <vector>

namespace index_converter::s_squared {
class Level {
  public:
    Level(
        std::shared_ptr<std::vector<spin_algebra::Multiplicity>> initialMultiplicities,
        size_t number_of_summations);
    void setMultiplicity(size_t number, spin_algebra::Multiplicity multiplicity);

    std::shared_ptr<std::vector<spin_algebra::Multiplicity>> getInitialMultiplicities() const;

    std::strong_ordering operator<=>(const Level& rhs) const;

    spin_algebra::Multiplicity getMultiplicity(size_t number) const;
    double getSpin(size_t number) const;
    size_t getSize() const;
    spin_algebra::Multiplicity total() const;

  private:
    std::shared_ptr<std::vector<spin_algebra::Multiplicity>> initialMultiplicities_;
    std::vector<spin_algebra::Multiplicity> intermediateMultiplicities_;
};
} // namespace index_converter::s_squared

#endif // SPINNER_S2LEVEL_H