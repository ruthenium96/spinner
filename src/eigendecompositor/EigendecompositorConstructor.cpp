#include "EigendecompositorConstructor.h"

#include "src/eigendecompositor/ExactEigendecompositor.h"
#include "src/eigendecompositor/ExplicitQuantitiesEigendecompositor.h"
#include "src/eigendecompositor/ImplicitSSquareEigendecompositor.h"
#include "src/eigendecompositor/OneSymbolInHamiltonianEigendecompositor.h"

namespace eigendecompositor {
std::unique_ptr<AbstractEigendecompositor> EigendecompositorConstructor::construct(
    const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
    const quantum::linear_algebra::FactoriesList& factories) {

    const auto& indexConverter =
        consistentModelOptimizationList.getModel().getIndexConverter();

    std::unique_ptr<eigendecompositor::AbstractEigendecompositor> eigendecompositor =
        std::make_unique<eigendecompositor::ExactEigendecompositor>(
            indexConverter,
            factories);

    const auto& symbolic_worker =
        consistentModelOptimizationList.getModel().getSymbolicWorker();

    // todo: we need only J and D if there is no field
    size_t number_of_all_J =
        symbolic_worker.getAllNames(model::symbols::J).size();
    size_t number_of_all_D =
        symbolic_worker.getAllNames(model::symbols::D).size();
    if (number_of_all_J == 1 && number_of_all_D == 0
        || number_of_all_J == 0 && number_of_all_D == 1) {
        model::symbols::SymbolName symbol_name;
        if (number_of_all_J == 1) {
            symbol_name = symbolic_worker.getChangeableNames(model::symbols::J)[0];
        } else if (number_of_all_D == 1) {
            symbol_name = symbolic_worker.getChangeableNames(model::symbols::D)[0];
        }
        auto getter = [&symbolic_worker, symbol_name]() {
            return symbolic_worker.getValueOfName(symbol_name);
        };
        eigendecompositor =
            std::make_unique<eigendecompositor::OneSymbolInHamiltonianEigendecompositor>(
                std::move(eigendecompositor),
                getter);
    }

    if (consistentModelOptimizationList.getOptimizationList().isSSquaredTransformed()) {
        eigendecompositor = std::make_unique<eigendecompositor::ImplicitSSquareEigendecompositor>(
            std::move(eigendecompositor),
            factories);
    }

    eigendecompositor = std::make_unique<eigendecompositor::ExplicitQuantitiesEigendecompositor>(
        std::move(eigendecompositor),
        indexConverter,
        factories);

    return std::move(eigendecompositor);
}

}  // namespace eigendecompositor