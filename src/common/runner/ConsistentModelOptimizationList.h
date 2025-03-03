#ifndef SPINNER_CONSISTENTMODELOPTIMIZATIONLIST_H
#define SPINNER_CONSISTENTMODELOPTIMIZATIONLIST_H

#include <memory>
#include "src/common/physical_optimization/OptimizationList.h"
#include "src/model/Model.h"

#include "src/common/index_converter/lexicographic/IndexConverter.h"
#include "src/common/index_converter/s_squared/IndexConverter.h"

namespace runner {
// This class keeps consistent pair (Model, OptimizationList).
// If inconsistent pair was passed to constructor, it throws a suitable exception.
class ConsistentModelOptimizationList {
  public:
    ConsistentModelOptimizationList(
        model::ModelInput,
        common::physical_optimization::OptimizationList);

    void InitializeDerivatives();
    void
    setNewValueToChangeableSymbol(const model::symbols::SymbolName& symbol_name, double new_value);

    const model::Model& getModel() const;
    const common::physical_optimization::OptimizationList& getOptimizationList() const;

    std::shared_ptr<const index_converter::lexicographic::IndexConverter> getLexIndexConverter() const;
    std::shared_ptr<const index_converter::AbstractIndexConverter> getIndexConverter() const;
    std::shared_ptr<const index_converter::s_squared::IndexConverter> getSquareIndexConveter() const;

    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
    getOperatorsForExplicitConstruction() const;
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>&
    getDerivativeOperatorsForExplicitConstruction() const;

    bool isImplicitSSquarePossible() const;
    bool isImplicitMSquarePossible() const;
    bool isExplicitMSquarePossible() const;
    bool isGSquaredT00Possible() const;
    bool isGSzSquaredPossible() const;

  private:
    std::unique_ptr<model::Model> model_;
    common::physical_optimization::OptimizationList optimizationList_;
    std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>
        operators_for_explicit_construction_;
    std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>
        derivatives_for_explicit_construction_;

    // TODO: replace it with abstract factory
    std::shared_ptr<index_converter::s_squared::IndexConverter> s_squared_index_converter_;
    std::shared_ptr<index_converter::lexicographic::IndexConverter> lex_index_converter_;
};
}  // namespace runner

#endif  //SPINNER_CONSISTENTMODELOPTIMIZATIONLIST_H
