#include "OptimizedSpaceConstructor.h"

#include "src/space/optimization/NonAbelianSimplifier.h"
#include "src/space/optimization/PositiveProjectionsEliminator.h"
#include "src/space/optimization/S2Transformer.h"
#include "src/space/optimization/Symmetrizer.h"
#include "src/space/optimization/TzSorter.h"
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

        //        if (space_history_.isNonAbelianSimplified && !new_group.properties.is_abelian) {
        //            throw std::invalid_argument(
        //                "Symmetrization after using of non-Abelian simplifier causes bugs.");
        //        }

        Symmetrizer symmetrizer(indexConverter, group, factories);
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
    //    space::optimization::NonAbelianSimplifier nonAbelianSimplifier(factories);
    //    space_ = nonAbelianSimplifier.apply(std::move(space_));
    //    space_history_.number_of_non_simplified_abelian_groups = 0;
    //    space_history_.isNonAbelianSimplified = true;

    if (optimizationList.isSSquaredTransformed()) {
        const auto number_of_mults = indexConverter.get_mults().size();
        // TODO: fix it somehow based on the information from groups:
        const size_t number_of_summation = number_of_mults - 1;

        std::vector<std::vector<std::set<size_t>>> all_groups_orbits_of_mults;
        for (const auto& group : optimizationList.getGroupsToApply()) {
            all_groups_orbits_of_mults.emplace_back(group.construct_orbits_of_mults());
        }

        // TODO: can we move it to constructor of S2Transformer?
        auto order_of_summations = spin_algebra::OrderOfSummation::constructFromOrbits(
            all_groups_orbits_of_mults,
            number_of_mults,
            number_of_summation);

        S2Transformer transformer(indexConverter, factories, order_of_summations);
        space = transformer.apply(std::move(space));
    }

    if (!spaceIsNormalized) {
        for (auto& subspace : space.getBlocks()) {
            subspace.decomposition->normalize();
        }
    }

    return space;
}
}  // namespace space::optimization