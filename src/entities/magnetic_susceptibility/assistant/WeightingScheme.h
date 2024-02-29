#ifndef SPINNER_WEIGHTINGSCHEME_H
#define SPINNER_WEIGHTINGSCHEME_H

#include <vector>

namespace magnetic_susceptibility {

enum WeightingSchemeEnum { per_point, per_interval };

class WeightingScheme {
  public:
    WeightingScheme() = default;
    WeightingScheme(std::vector<double> abscissa, WeightingSchemeEnum weightingSchemeEnum);
    double at(size_t i) const;
    const std::vector<double>& getWeights() const;

  private:
    std::vector<double> weights_;
};

}  // namespace magnetic_susceptibility

#endif  //SPINNER_WEIGHTINGSCHEME_H
