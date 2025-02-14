#include "OptimizedSpaceConstructor.h"
#include <memory>

#include "src/common/Logger.h"
#include "src/common/physical_optimization/OptimizationList.h"
#include "src/space/optimization/NonAbelianSimplifier.h"
#include "src/space/optimization/PositiveProjectionsEliminator.h"
#include "src/space/optimization/NonMinimalProjectionsEliminator.h"
#include "src/space/optimization/S2Transformer.h"
#include "src/space/optimization/Symmetrizer.h"
#include "src/space/optimization/TSquaredSorter.h"
#include "src/space/optimization/TzSorter.h"
#include "src/spin_algebra/GroupAdapter.h"
#include "src/common/index_converter/s_squared/OrderOfSummation.h"
#include "src/common/index_converter/lexicographic/IndexPermutator.h"
#include "src/common/index_converter/s_squared/IndexPermutator.h"


namespace {
std::vector<size_t> sizes_of_blocks(const space::Space& space) {
    std::vector<size_t> answer(space.getBlocks().size());
    for (int i = 0; i < space.getBlocks().size(); ++i) {
        answer[i] = space.getBlocks()[i].decomposition->size_cols();
    }
    return answer;
}
}

namespace space::optimization {

// TODO: if Space/Subspace will be transformed into classes,
//  it should be static constructor of Space.
Space OptimizedSpaceConstructor::construct(
    const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
    const quantum::linear_algebra::FactoriesList& factories) {
    std::shared_ptr<const lexicographic::IndexConverter> indexConverter =
        consistentModelOptimizationList.getModel().getIndexConverter();
    const common::physical_optimization::OptimizationList& optimizationList =
        consistentModelOptimizationList.getOptimizationList();

    Space space = Space(
        indexConverter->get_total_space_size(),
        factories);

    bool spaceIsNormalized = true;

    if (optimizationList.isTzSorted()) {
        TzSorter tz_sorter(indexConverter, factories);
        common::Logger::detailed_msg("Tz-sortation has started.");
        space = tz_sorter.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("Tz-sortation is finished.");
    }

    common::Logger::separate(1, common::PrintLevel::detailed);

    if (optimizationList.isITOBasis() && optimizationList.isTSquaredSorted()) {
        TSquaredSorter tsquared_sorter(sSquaredIndexConverter, factories);
        common::Logger::detailed_msg("T2-sortation has started.");
        space = tsquared_sorter.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("T2-sortation is finished.");
    }

    if (optimizationList.isPositiveProjectionsEliminated()) {
        uint32_t max_ntz_proj = indexConverter->get_max_ntz_proj();

        PositiveProjectionsEliminator positiveProjectionsEliminator(max_ntz_proj);
        common::Logger::detailed_msg("Positive projections elimination has started.");
        space = positiveProjectionsEliminator.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("Positive projections elimination is finished.");
        common::Logger::separate(1, common::detailed);
    }

    if (optimizationList.isNonMinimalProjectionsEliminated() && optimizationList.isITOBasis()) {
        uint32_t max_ntz_proj = indexConverter->get_max_ntz_proj();

        NonMinimalProjectionsEliminator nonMinimalProjectionsEliminator(max_ntz_proj);
        common::Logger::detailed_msg("Non-minimal projections elimination has started.");
        space = nonMinimalProjectionsEliminator.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("Non-minimal projections elimination is finished.");
        common::Logger::separate(1, common::detailed);
    }

    for (size_t i = 0; i < optimizationList.getGroupsToApply().size(); ++i) {
        // Symmetrization breaks normalization, because of using integer values of coefficients.
        spaceIsNormalized = false;

        const auto& group = optimizationList.getGroupsToApply().at(i);
        std::shared_ptr<const index_converter::AbstractIndexPermutator> permutator;
        if (optimizationList.isITOBasis()) {
            permutator = 
                std::make_shared<index_converter::s_squared::IndexPermutator>(sSquaredIndexConverter, group);
        } else {
            permutator = 
                std::make_shared<index_converter::lexicographic::IndexPermutator>(lexIndexConverter, group);
        }

        Symmetrizer symmetrizer(permutator, group, factories);
        common::Logger::detailed_msg("Symmetrization has started.");
        space = symmetrizer.apply(std::move(space));
        common::Logger::verbose("Sizes of blocks:\n{}", fmt::join(sizes_of_blocks(space), ", "));
        common::Logger::detailed_msg("Symmetrization is finished.");
        if (i + 1 == optimizationList.getGroupsToApply().size()) {
            common::Logger::separate(1, common::detailed);
        } else {
            common::Logger::separate(2, common::verbose);
        }
    }

    if (!spaceIsNormalized) {
        for (auto& subspace : space.getBlocks()) {
            subspace.decomposition->normalize();
        }
    }

    if (optimizationList.isSSquaredTransformed()) {

        S2Transformer transformer(
            lexIndexConverter,
            factories,
            consistentModelOptimizationList.getSSquaredConverter());

        common::Logger::detailed_msg("S2-transformation has started.");
        common::orderOfSummationPrint(
            *consistentModelOptimizationList.getSSquaredConverter()->getOrderOfSummation());
        space = transformer.apply(std::move(space));
        common::Logger::detailed_msg("S2-transformation is finished.");
    }

    common::Logger::separate(0, common::detailed);

    return space;
}
}  // namespace space::optimization