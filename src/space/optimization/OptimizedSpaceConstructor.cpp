#include "OptimizedSpaceConstructor.h"

#include "src/space/optimization/NonAbelianSimplifier.h"
#include "src/space/optimization/PositiveProjectionsEliminator.h"
#include "src/space/optimization/Symmetrizer.h"
#include "src/space/optimization/TzSorter.h"

namespace space::optimization {

Space OptimizedSpaceConstructor::construct(
    const runner::ConsistentModelOptimizationList& consistentModelOptimizationList) {
    const lexicographic::IndexConverter& indexConverter =
        consistentModelOptimizationList.getModel().getIndexConverter();
    const common::physical_optimization::OptimizationList& optimizationList =
        consistentModelOptimizationList.getOptimizationList();

    Space space = Space(
        consistentModelOptimizationList.getModel().getIndexConverter().get_total_space_size());

    bool spaceIsNormalized = true;

    if (optimizationList.isTzSorted()) {
        TzSorter tz_sorter(indexConverter);
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

        //        if (space_history_.isNonAbelianSimplified && !new_group.properties.is_abelian) {
        //            throw std::invalid_argument(
        //                "Symmetrization after using of non-Abelian simplifier causes bugs.");
        //        }

        space::optimization::Symmetrizer symmetrizer(indexConverter, group);
        space = symmetrizer.apply(std::move(space));

        //        if (!new_group.properties.is_abelian) {
        //            ++space_history_.number_of_non_simplified_abelian_groups;
        //        }
    }
    //    if (space_history_.number_of_non_simplified_abelian_groups == 0) {
    //        return;
    //    }
    //    if (space_history_.number_of_non_simplified_abelian_groups != 1) {
    //        throw std::invalid_argument(
    //            "Non-Abelian simplification after using of two Non-Abelian Symmetrizers "
    //            "currently is not allowed.");
    //    }
    //    space::optimization::NonAbelianSimplifier nonAbelianSimplifier;
    //    space_ = nonAbelianSimplifier.apply(std::move(space_));
    //    space_history_.number_of_non_simplified_abelian_groups = 0;
    //    space_history_.isNonAbelianSimplified = true;

    if (!spaceIsNormalized) {
        for (auto& subspace : space.getBlocks()) {
            // TODO: maybe, we can implement normalize as Space method
            subspace.decomposition->normalize();
        }
    }

    return space;
}
}  // namespace space::optimization