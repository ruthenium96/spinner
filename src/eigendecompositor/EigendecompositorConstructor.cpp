#include "EigendecompositorConstructor.h"

#include "src/common/Logger.h"
#include "src/common/Quantity.h"
#include "src/eigendecompositor/ExactEigendecompositor.h"
#include "src/eigendecompositor/ExplicitQuantitiesEigendecompositor.h"
#include "src/eigendecompositor/FTLMEigendecompositor.h"
#include "src/eigendecompositor/ImplicitQuantityEigendecompositor.h"
#include "src/eigendecompositor/OneSymbolInHamiltonianEigendecompositor.h"

namespace eigendecompositor {
std::unique_ptr<AbstractEigendecompositor> EigendecompositorConstructor::construct(
    const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
    const quantum::linear_algebra::FactoriesList& factories) {

    const auto& indexConverter =
        consistentModelOptimizationList.getIndexConverter();

    common::Logger::detailed_msg("Eigendecompositor information:");
    common::Logger::detailed("ExactEigendecompositor will be used");
    std::unique_ptr<eigendecompositor::AbstractEigendecompositor> eigendecompositor;

    if (consistentModelOptimizationList.getOptimizationList().isFTLMApproximated()) {
        auto FTLM_settings = consistentModelOptimizationList.getOptimizationList().getFTLMSettings();
        common::Logger::detailed("FTLMEigendecompositor will be used");
        eigendecompositor =
            std::make_unique<eigendecompositor::FTLMEigendecompositor>(
                indexConverter,
                factories,
                FTLM_settings.krylov_subspace_size,
                FTLM_settings.exact_decomposition_threshold,
                FTLM_settings.number_of_seeds);
    } else {
        eigendecompositor = std::make_unique<eigendecompositor::ExactEigendecompositor>(
            indexConverter,
            factories);
    }

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
        common::Logger::detailed(
            "OneSymbolInHamiltonianEigendecompositor will be used, "
            "the name of the parameter is {}.",
            symbol_name.get_name());
        eigendecompositor =
            std::make_unique<eigendecompositor::OneSymbolInHamiltonianEigendecompositor>(
                std::move(eigendecompositor),
                getter);
    }

    if (consistentModelOptimizationList.isImplicitSSquarePossible()) {
        common::Logger::detailed("ImplicitSSquareEigendecompositor will be used for S_total_squared.");
        eigendecompositor = std::make_unique<eigendecompositor::ImplicitQuantityEigendecompositor>(
            std::move(eigendecompositor),
            factories,
            common::S_total_squared,
            indexConverter->get_max_ntz_proj());    
    }
    if (consistentModelOptimizationList.isImplicitMSquarePossible()) {
            common::Logger::detailed("ImplicitSSquareEigendecompositor will be used for M_total_squared.");
            eigendecompositor = std::make_unique<eigendecompositor::ImplicitQuantityEigendecompositor>(
                std::move(eigendecompositor),
                factories,
                common::M_total_squared,
                indexConverter->get_max_ntz_proj());
    }

    common::Logger::detailed("ExplicitQuantitiesEigendecompositor will be used.");
    eigendecompositor = std::make_unique<eigendecompositor::ExplicitQuantitiesEigendecompositor>(
        std::move(eigendecompositor),
        indexConverter,
        factories);

    common::Logger::separate(0, common::detailed);

    return std::move(eigendecompositor);
}

}  // namespace eigendecompositor