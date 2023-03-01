#ifndef SPINNER_MODELOPTIMIZATIONLISTCONSISTENCE_H
#define SPINNER_MODELOPTIMIZATIONLISTCONSISTENCE_H

#include "src/common/physical_optimization/OptimizationList.h"
#include "src/model/Model.h"

namespace runner {
// This class throws exception if Model and OptimizationList are inconsistent.
class ModelOptimizationListConsistence {
  public:
    static void check(
        const model::Model& model,
        const common::physical_optimization::OptimizationList& optimizationList);

  private:
    // This function checks if multiplicities are invariant to all elements from the group.
    static void
    checkMultiplicitiesGroupConsistence(const std::vector<int>&, const group::Group& group);
    // This function checks if all SymbolNames from SymbolicWorker are invariant to all elements from the group.
    static void checkAllSymbolNamesGroupConsistence(
        const model::symbols::SymbolicWorker& symbols,
        const group::Group& group);
    // This function checks if two-center SymbolNames set are invariant to element permutation.
    // It requires getter from SymbolicWorker, for example, getIsotropicExchangeSymbolName.
    static bool checkSymbolNamesGroupElementConsistence(
        const std::function<std::optional<model::symbols::SymbolName>(size_t, size_t j)>& getter,
        group::Permutation element);
    // This function checks if one-center SymbolNames set are invariant to element permutation.
    // It requires getter from SymbolicWorker, for example, getGFactorSymbolName.
    static bool checkSymbolNamesGroupElementConsistence(
        const std::function<model::symbols::SymbolName(size_t)>& getter,
        group::Permutation element);
};
}  // namespace runner
#endif  //SPINNER_MODELOPTIMIZATIONLISTCONSISTENCE_H
