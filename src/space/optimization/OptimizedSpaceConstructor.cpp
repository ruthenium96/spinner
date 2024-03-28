#include "OptimizedSpaceConstructor.h"

#include "src/space/optimization/NonAbelianSimplifier.h"
#include "src/space/optimization/PositiveProjectionsEliminator.h"
#include "src/space/optimization/S2Transformer.h"
#include "src/space/optimization/Symmetrizer.h"
#include "src/space/optimization/TzSorter.h"
#include "src/spin_algebra/GroupAdapter.h"
#include "src/spin_algebra/OrderOfSummation.h"

namespace space::optimization {

// TODO: if Space/Subspace will be transformed into classes,
//  it should be static constructor of Space.
Space OptimizedSpaceConstructor::construct(
    const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
    const quantum::linear_algebra::FactoriesList& factories) {
    const lexicographic::IndexConverter& indexConverter =
        consistentModelOptimizationList.getModel().getIndexConverter();
    const common::physical_optimization::OptimizationList& optimizationList =
        consistentModelOptimizationList.getOptimizationList();

    Space space = Space(
        consistentModelOptimizationList.getModel().getIndexConverter().get_total_space_size(),
        factories);

    bool spaceIsNormalized = true;

    if (optimizationList.isTzSorted()) {
        TzSorter tz_sorter(indexConverter, factories);
        space = tz_sorter.apply(std::move(space));
    }
    if (optimizationList.isPositiveProjectionsEliminated()) {
        uint32_t max_ntz_proj = indexConverter.get_max_ntz_proj();

        PositiveProjectionsEliminator positiveProjectionsEliminator(max_ntz_proj);
        space = positiveProjectionsEliminator.apply(std::move(space));
    }
    for (const auto& group : optimizationList.getGroupsToApply()) {
        // Symmetrization breaks normalization, because of using integer values of coefficients.
        spaceIsNormalized = false;

        Symmetrizer symmetrizer(indexConverter, group, factories);
        space = symmetrizer.apply(std::move(space));

        //        if (!new_group.properties.is_abelian) {
        //            ++space_history_.number_of_non_simplified_abelian_groups;
        //        }
    }

    if (!spaceIsNormalized) {
        for (auto& subspace : space.getBlocks()) {
            subspace.decomposition->normalize();
        }
    }

    if (optimizationList.isSSquaredTransformed()) {
        const auto number_of_mults = indexConverter.get_mults().size();
        auto group_adapter =
            spin_algebra::GroupAdapter(optimizationList.getGroupsToApply(), number_of_mults);

        S2Transformer transformer(
            indexConverter,
            factories,
            group_adapter.getOrderOfSummations(),
            group_adapter.getRepresentationMultiplier());
        space = transformer.apply(std::move(space));
    }

    return space;
}
}  // namespace space::optimization