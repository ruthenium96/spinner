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
    initialize();
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        const auto& subspace = space.getBlocks().at(i);
        auto operators_to_calculate = operators;
        auto derivatives_operators_to_calculate = derivatives_operators;
        BuildSubspectra(operators_to_calculate, derivatives_operators_to_calculate, i, subspace);
        assert(operators_to_calculate.size() == 0);
        assert(derivatives_operators_to_calculate.size() == 0);
    }
    finalize();
}
}  // namespace eigendecompositor