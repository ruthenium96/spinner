#ifndef SPINNER_TSQUAREDSORTER_H
#define SPINNER_TSQUAREDSORTER_H

#include "src/common/index_converter/s_squared/IndexConverter.h"
#include "src/entities/data_structures/FactoriesList.h"
#include "src/space/Space.h"
#include "src/spin_algebra/Multiplicity.h"

namespace space::optimization {
class TSquaredSorter {
public:
    TSquaredSorter(std::shared_ptr<const index_converter::s_squared::IndexConverter> indexConverter,
        quantum::linear_algebra::FactoriesList factories);

    Space apply(Space&& space) const;

private:
    std::shared_ptr<const index_converter::s_squared::IndexConverter> converter_;
    const quantum::linear_algebra::FactoriesList factories_;
    uint32_t max_total_mult_;

    spin_algebra::Multiplicity number_of_block_to_multiplicity(size_t number_of_block) const;
    size_t multiplicity_to_number_of_block(spin_algebra::Multiplicity multiplicity) const;
};
}  // namespace space::optimization

#endif // SPINNER_TSQUAREDSORTER_H