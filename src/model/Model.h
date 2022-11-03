#ifndef SPINNER_MODEL_H
#define SPINNER_MODEL_H

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
    Model& InitializeZeroFieldSplitting();
    Model& InitializeSSquared();
    Model& InitializeGSzSquared();
    Model& InitializeGSzSquaredDerivatives();

    bool is_s_squared_initialized() const;
    bool is_g_sz_squared_initialized() const;
    bool is_g_sz_squared_derivatives_initialized() const;
    bool is_isotropic_exchange_derivatives_initialized() const;
    bool is_zero_field_splitting_initialized() const;

    const lexicographic::IndexConverter& getIndexConverter() const;
    const symbols::Symbols& getSymbols() const;
    symbols::Symbols& getSymbols();
    const operators::Operator& getOperator(common::QuantityEnum) const;
    const operators::Operator& getOperatorDerivative(
        common::QuantityEnum,
        const symbols::SymbolName&) const;

  private:
    struct OperatorsHistory {
        bool isotropic_exchange_in_hamiltonian = false;
        bool s_squared = false;
        bool g_sz_squared = false;
        bool isotropic_exchange_derivatives = false;
        bool g_sz_squared_derivatives = false;
        bool zfs_in_hamiltonian = false;
    };
    // TODO: should we move it from here?
    const lexicographic::IndexConverter converter_;
    symbols::Symbols symbols_;

    // std::shared_ptr<BasicQuantity>, where BasicQuantity is virtual class?
    // or Quantity can consist std::unique_ptr<VirtualMatrix>?
    // or can we just separate Operator and Spectrum?
    operators::Operator energy_operator;
    std::optional<operators::Operator> s_squared_operator;
    std::optional<operators::Operator> g_sz_squared_operator;
    std::map<std::pair<common::QuantityEnum, symbols::SymbolName>, operators::Operator>
        derivatives_map_;

    OperatorsHistory operators_history_;
};
}  // namespace model

#endif  //SPINNER_MODEL_H
