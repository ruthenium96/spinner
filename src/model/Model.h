#ifndef JULY_MODEL_H
#define JULY_MODEL_H

#include <vector>

#include "src/common/Quantity.h"
#include "src/model/operators/Operator.h"
#include "src/model/symbols/Symbols.h"
namespace model {
// class Model is responsible for Operators (parametrized by Symbols) and multiplicities
class Model {
  public:
    explicit Model(const std::vector<int>& mults);

    Model& InitializeIsotropicExchange();
    Model& InitializeIsotropicExchangeDerivatives();
    Model& InitializeSSquared();

    bool is_s_squared_initialized() const;
    bool is_isotropic_exchange_derivatives_initialized() const;

    const lexicographic::IndexConverter& getIndexConverter() const;
    const symbols::Symbols& getSymbols() const;
    symbols::Symbols& getSymbols();
    const operators::Operator& getOperator(common::QuantityEnum) const;
    const operators::Operator& getOperatorDerivative(
        common::QuantityEnum,
        symbols::SymbolTypeEnum,
        const symbols::SymbolName&) const;

  private:
    struct OperatorsHistory {
        bool isotropic_exchange_in_hamiltonian = false;
        bool s_squared = false;
        bool isotropic_exchange_derivatives = false;
    };
    const lexicographic::IndexConverter converter_;
    symbols::Symbols symbols_;

    // std::shared_ptr<BasicQuantity>, where BasicQuantity is virtual class?
    // or Quantity can consist std::unique_ptr<VirtualMatrix>?
    // or can we just separate Operator and Spectrum?
    operators::Operator energy_operator;
    std::optional<operators::Operator> s_squared_operator;
    std::map<symbols::SymbolName, operators::Operator>
        derivative_of_energy_wrt_exchange_parameters_operator;

    OperatorsHistory operators_history_;
};
}  // namespace model

#endif  //JULY_MODEL_H
