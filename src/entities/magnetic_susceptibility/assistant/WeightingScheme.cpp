#include "WeightingScheme.h"

#include <numeric>

#include "src/common/IncorrectEnumError.h"

namespace magnetic_susceptibility {

WeightingScheme::WeightingScheme(
    std::vector<double> abscissa,
    WeightingSchemeEnum weightingSchemeEnum) {
    if (weightingSchemeEnum == per_point) {
        weights_.resize(abscissa.size(), 1);
    } else if (weightingSchemeEnum == per_interval) {
        if (abscissa.size() == 1) {
            throw std::logic_error(
                "Cannot apply per_interval weighting scheme for only one point.");
        }
        weights_.resize(abscissa.size());
        // for all points, except the leftest and the rightest...
        for (size_t i = 1; i < abscissa.size() - 1; ++i) {
            double diff_one = abscissa[i] - abscissa[i - 1];
            double diff_two = abscissa[i + 1] - abscissa[i];
            weights_[i] = (diff_one + diff_two) / 2;
        }
        // the leftest point:
        double diff_one = abscissa[1] - abscissa[0];
        weights_[0] = diff_one;
        // the rightest point:
        double diff_two = abscissa[abscissa.size() - 1] - abscissa[abscissa.size() - 2];
        weights_[abscissa.size() - 1] = diff_two;

        double total_weight = std::accumulate(weights_.begin(), weights_.end(), 0.0);
        for (auto& weight : weights_) {
            weight /= total_weight;
            // reweight to make similar to per_point
            weight *= abscissa.size();
        }
    } else {
        throw common::IncorrectEnumError(
            "Unknown type of WeightingSchemeEnum has been passed to WeightingScheme");
    }
}

double WeightingScheme::at(size_t i) const {
    return weights_.at(i);
}

const std::vector<double>& WeightingScheme::getWeights() const {
    return weights_;
}
}  // namespace magnetic_susceptibility