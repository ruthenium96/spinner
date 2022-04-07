#include "ModelOptimizationListConsistence.h"

namespace runner {
void ModelOptimizationListConsistence::check(
    const model::Model& model,
    const common::physical_optimization::OptimizationList& optimizationList) {
    // Symmetrizer can be applied for all Hamiltonian terms,
    // if they match up with symmetry:
    for (const auto& group : optimizationList.getGroupsToApply()) {
        checkAllSymbolNamesGroupConsistence(model.getSymbols(), group);
    }
}

void ModelOptimizationListConsistence::checkAllSymbolNamesGroupConsistence(
    const model::symbols::Symbols& symbols,
    const group::Group& group) {
    if (symbols.isIsotropicExchangeInitialized()) {
        for (const auto& element : group.getElements()) {
            std::function<model::symbols::SymbolName(size_t, size_t)> getterFunction =
                [ObjectPtr = &symbols](size_t i, size_t j) {
                    return ObjectPtr->getIsotropicExchangeSymbolName(i, j);
                };
            if (!checkSymbolNamesGroupElementConsistence(getterFunction, element)) {
                throw std::invalid_argument(
                    "Isotropic exchange symbols do not match applied symmetries");
                // TODO: should we rename this exception?
            }
        }
    }
    if (symbols.isGFactorInitialized()) {
        for (const auto& element : group.getElements()) {
            std::function<model::symbols::SymbolName(size_t)> getterFunction =
                [ObjectPtr = &symbols](size_t i) { return ObjectPtr->getGFactorSymbolName(i); };
            if (!checkSymbolNamesGroupElementConsistence(getterFunction, element)) {
                throw std::invalid_argument("g factor symbols do not match applied symmetries");
                // TODO: should we rename this exception?
            }
        }
    }
}

bool ModelOptimizationListConsistence::checkSymbolNamesGroupElementConsistence(
    const std::function<model::symbols::SymbolName(size_t, size_t)>& getter,
    group::Permutation element) {
    // Construct here both initialSymbols and permutatedSymbols:
    std::vector<std::vector<model::symbols::SymbolName>> initialSymbols(
        element.size(),
        std::vector<model::symbols::SymbolName>(element.size()));
    auto permutatedSymbols = initialSymbols;
    for (size_t i = 0; i < element.size(); ++i) {
        for (size_t j = 0; j < element.size(); ++j) {
            auto symbolName = getter(i, j);
            initialSymbols[i][j] = symbolName;
            permutatedSymbols[element[i]][element[j]] = symbolName;
        }
    }
    return initialSymbols == permutatedSymbols;
}

bool ModelOptimizationListConsistence::checkSymbolNamesGroupElementConsistence(
    const std::function<model::symbols::SymbolName(size_t)>& getter,
    group::Permutation element) {
    std::vector<model::symbols::SymbolName> initialSymbols(element.size());
    auto permutatedSymbols = initialSymbols;

    for (size_t i = 0; i < element.size(); ++i) {
        auto symbolName = getter(i);
        initialSymbols[i] = symbolName;
        permutatedSymbols[element[i]] = symbolName;
    }
    return initialSymbols == permutatedSymbols;
}
}  // namespace runner
