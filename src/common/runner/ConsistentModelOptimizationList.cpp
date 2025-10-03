#include "ConsistentModelOptimizationList.h"

#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>

#include "src/common/Quantity.h"
#include "src/common/index_converter/lexicographic/IndexConverter.h"
#include "src/common/index_converter/s_squared/IndexConverter.h"
#include "src/common/physical_optimization/OptimizationList.h"
#include "src/model/operators/AbstractOperatorConstructor.h"
#include "src/model/operators/ITOOperatorConstructor.h"
#include "src/model/operators/LexOperatorConstructor.h"
#include "src/spin_algebra/GroupAdapter.h"

namespace {
void checkMultiplicitiesGroupConsistence(
    const std::vector<spin_algebra::Multiplicity>& mults,
    const group::Group& group);
bool checkSymbolNamesGroupElementConsistence(
    const std::function<std::optional<model::symbols::SymbolName>(size_t, size_t)>& getter,
    group::Permutation element);
bool checkSymbolNamesGroupElementConsistence(
    const std::function<model::symbols::SymbolName(size_t)>& getter,
    group::Permutation element);
void checkAllSymbolNamesGroupConsistence(
    const model::symbols::SymbolicWorker& symbols,
    const group::Group& group);
void checkModelOptimizationListConsistence(
    const model::ModelInput& modelInput,
    const common::physical_optimization::OptimizationList& optimizationList);
}  // namespace

namespace runner {

const model::Model& ConsistentModelOptimizationList::getModel() const {
    return *model_;
}

const common::physical_optimization::OptimizationList&
ConsistentModelOptimizationList::getOptimizationList() const {
    return optimizationList_;
}

ConsistentModelOptimizationList::ConsistentModelOptimizationList(
    model::ModelInput modelInput,
    common::physical_optimization::OptimizationList optimizationList) :
    optimizationList_(std::move(optimizationList)) {
    // this function throw an exception if ModelInput and OptimizationList are inconsistent
    checkModelOptimizationListConsistence(modelInput, optimizationList_);

    std::unique_ptr<model::operators::AbstractOperatorConstructor> operator_constructor;
    if (getOptimizationList().isITOBasis()) {
        const auto number_of_mults = modelInput.getMults().size();
        auto group_adapter =
            spin_algebra::GroupAdapter(optimizationList_.getGroupsToApply(), number_of_mults);

        s_squared_index_converter_ = std::make_unique<index_converter::s_squared::IndexConverter>(
            modelInput.getMults(), 
            group_adapter.getOrderOfSummations());

        operator_constructor = std::make_unique<model::operators::ITOOperatorConstructor>(s_squared_index_converter_);
    } else {
        lex_index_converter_ = std::make_shared<index_converter::lexicographic::IndexConverter>(modelInput.getMults());
        operator_constructor = std::make_unique<model::operators::LexOperatorConstructor>(lex_index_converter_);
    }

    model_ = std::make_unique<model::Model>(std::move(modelInput), std::move(operator_constructor));

    operators_for_explicit_construction_[common::Energy] =
        model_->getOperator(common::Energy).value();

    if (isImplicitSSquarePossible() || isImplicitMSquarePossible()) {
        return;
    }
    if (isExplicitMSquarePossible()) {
        operators_for_explicit_construction_[common::M_total_squared] = 
            model_->getOperator(common::M_total_squared).value();
        return;
    }
    if (isGSquaredT00Possible()) {
        operators_for_explicit_construction_[common::g_squared_T00] = 
            model_->getOperator(common::g_squared_T00).value();
        return;
    }
    if (isGSzSquaredPossible()) {
        operators_for_explicit_construction_[common::gSz_total_squared] =
            model_->getOperator(common::gSz_total_squared).value();
        return;
    }
}

void ConsistentModelOptimizationList::InitializeDerivatives() {
    if (derivatives_for_explicit_construction_.empty()) {
        for (const auto& [pair, shared_ptr] : model_->getOperatorDerivatives()) {
            auto quantity_enum = pair.first;
            auto type_of_symbol = getModel().getSymbolicWorker().getSymbolProperty(pair.second).type_enum.value();
            if (type_of_symbol == model::symbols::g_factor) {
                if (isImplicitSSquarePossible() || isImplicitMSquarePossible()) {
                    continue;
                }
                if (isExplicitMSquarePossible() && quantity_enum != common::M_total_squared) {
                    continue;
                }
                if (isGSquaredT00Possible() && quantity_enum != common::g_squared_T00) {
                    continue;
                }
                if (isGSzSquaredPossible() && quantity_enum != common::gSz_total_squared) {
                    continue;
                }
            }
            derivatives_for_explicit_construction_[pair] = shared_ptr;
        }
    }
}

void ConsistentModelOptimizationList::setNewValueToChangeableSymbol(
    const model::symbols::SymbolName& symbol_name,
    double new_value) {
    model_->setNewValueToChangeableSymbol(symbol_name, new_value);
}

const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
ConsistentModelOptimizationList::getOperatorsForExplicitConstruction() const {
    return operators_for_explicit_construction_;
}

const std::map<
    std::pair<common::QuantityEnum, model::symbols::SymbolName>,
    std::shared_ptr<const model::operators::Operator>>&
ConsistentModelOptimizationList::getDerivativeOperatorsForExplicitConstruction() const {
    return derivatives_for_explicit_construction_;
}

std::shared_ptr<const index_converter::lexicographic::IndexConverter>
ConsistentModelOptimizationList::getLexIndexConverter() const {
    return lex_index_converter_;
}

std::shared_ptr<const index_converter::s_squared::IndexConverter> 
ConsistentModelOptimizationList::getSquareIndexConveter() const {
    return s_squared_index_converter_;
}

std::shared_ptr<const index_converter::AbstractIndexConverter>
ConsistentModelOptimizationList::getIndexConverter() const {
    if (getSquareIndexConveter() != nullptr) {
        return getSquareIndexConveter();
    } else {
        return getLexIndexConverter();
    }
}

bool ConsistentModelOptimizationList::isImplicitSSquarePossible() const {
    return getModel().gFactorsAreAllNoneOrAreTheSame()
    && !getModel().getSymbolicWorker().isZFSInitialized()
    && getOptimizationList().isNonMinimalProjectionsEliminated();
}

bool ConsistentModelOptimizationList::isImplicitMSquarePossible() const {
    return getModel().gFactorsAreAllNoneOrAreTheSame()
    && !getOptimizationList().isNonMinimalProjectionsEliminated()
    && getOptimizationList().isTzSorted();
}

bool ConsistentModelOptimizationList::isExplicitMSquarePossible() const {
    return getModel().gFactorsAreAllNoneOrAreTheSame();
}

bool ConsistentModelOptimizationList::isGSquaredT00Possible() const {
    return !getModel().gFactorsAreNone() 
    && !getModel().getSymbolicWorker().isZFSInitialized()
    && getOptimizationList().isNonMinimalProjectionsEliminated();
}

bool ConsistentModelOptimizationList::isGSzSquaredPossible() const {
    return !getModel().gFactorsAreNone();
}

}  // namespace runner

namespace {
void checkModelOptimizationListConsistence(
    const model::ModelInput& modelInput,
    const common::physical_optimization::OptimizationList& optimizationList) {
    // Symmetrizer can be applied for all types of Hamiltonian terms...
    for (const auto& group : optimizationList.getGroupsToApply()) {
        // ...if spins invariant to group elements:
        checkMultiplicitiesGroupConsistence(modelInput.getMults(), group);
        // ...if Hamiltonian terms invariant to group elements:
        checkAllSymbolNamesGroupConsistence(modelInput.getSymbolicWorker(), group);
    }
    if (optimizationList.isSSquaredTransformed()) {
        // S2-transformation can be applied only for HDvV-Hamiltionian:
        if (modelInput.getSymbolicWorker().isZFSInitialized()) {
            throw std::invalid_argument("S2-transformation cannot be applied to ZFS-Hamiltonian");
        }
    }
    if (optimizationList.isLexBasis() && optimizationList.isTSquaredSorted() && !optimizationList.isSSquaredTransformed()) {
        throw std::invalid_argument("In the case of lex basis, we cannot TSquared Sort without SSquared Transformer");
    }
}

void checkMultiplicitiesGroupConsistence(
    const std::vector<spin_algebra::Multiplicity>& mults,
    const group::Group& group) {
    // TODO: split code above (maybe rewrite checkSymbolNamesGroupElementConsistence with templates?)
    if (mults.size() != group.size_of_permutations()) {
        throw std::length_error(
            "The size of group elements does not equal to the number of spins.");
    }
    for (const auto& el : group.getElements()) {
        std::vector<spin_algebra::Multiplicity> permutated_mults(mults);
        for (size_t i = 0; i < group.size_of_permutations(); ++i) {
            permutated_mults[i] = mults[el[i]];
        }
        if (permutated_mults != mults) {
            throw std::invalid_argument("Group permutes centers with different multiplicities.");
        }
    }
}

bool checkSymbolNamesGroupElementConsistence(
    const std::function<std::optional<model::symbols::SymbolName>(size_t, size_t)>& getter,
    group::Permutation element) {
    // Construct here both initialSymbols and permutatedSymbols:
    std::vector<std::vector<std::optional<model::symbols::SymbolName>>> initialSymbols(
        element.size(),
        std::vector<std::optional<model::symbols::SymbolName>>(element.size(), std::nullopt));
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

bool checkSymbolNamesGroupElementConsistence(
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

void checkAllSymbolNamesGroupConsistence(
    const model::symbols::SymbolicWorker& symbols,
    const group::Group& group) {
    if (symbols.isIsotropicExchangeInitialized()) {
        for (const auto& element : group.getElements()) {
            std::function<std::optional<model::symbols::SymbolName>(size_t, size_t)>
                getterFunction = [ObjectPtr = &symbols](size_t i, size_t j) {
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
}  // namespace