#include "AbstractEigendecompositor.h"

#include <cassert>

namespace eigendecompositor {

void AbstractEigendecompositor::BuildSpectra(
    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators,
    const space::Space& space) {
    buildSpectraWasCalled = true;
    {
        auto operators_to_calculate = operators;
        auto derivatives_operators_to_calculate = derivatives_operators;
        initialize(
            operators_to_calculate,
            derivatives_operators_to_calculate,
            space.getBlocks().size());
        assert(operators_to_calculate.size() == 0);
        assert(derivatives_operators_to_calculate.size() == 0);
    }

#pragma omp parallel for shared(space, operators, derivatives_operators) default(shared)
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        const auto& subspace = space.getBlocks().at(i);
        BuildSubspectra(i, subspace);
    }

    finalize();
}

bool AbstractEigendecompositor::BuildSpectraWasCalled() const {
    return buildSpectraWasCalled;
}
}  // namespace eigendecompositor