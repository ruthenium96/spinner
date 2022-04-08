#ifndef JULY_MODELOPTIMIZATIONLISTCONSISTENCE_H
#define JULY_MODELOPTIMIZATIONLISTCONSISTENCE_H

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
    // This function checks if all SymbolNames from Symbols are invariant to all elements from the group.
    static void checkAllSymbolNamesGroupConsistence(
        const model::symbols::Symbols& symbols,
        const group::Group& group);
    // This function checks if two-center SymbolNames set are invariant to element permutation.
    // It requires getter from Symbols, for example, getIsotropicExchangeSymbolName.
    static bool checkSymbolNamesGroupElementConsistence(
        const std::function<model::symbols::SymbolName(size_t, size_t j)>& getter,
        group::Permutation element);
    // This function checks if one-center SymbolNames set are invariant to element permutation.
    // It requires getter from Symbols, for example, getGFactorSymbolName.
    static bool checkSymbolNamesGroupElementConsistence(
        const std::function<model::symbols::SymbolName(size_t)>& getter,
        group::Permutation element);
};
}  // namespace runner
#endif  //JULY_MODELOPTIMIZATIONLISTCONSISTENCE_H
