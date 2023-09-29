#ifndef SPINNER_MODEL_H
#define SPINNER_MODEL_H

#include "ModelInput.h"
#include "symbols/NumericalWorker.h"

namespace model {

class Model {
  public:
    explicit Model(ModelInput modelInput);

    void InitializeDerivatives();

    // TODO: remove it?
    const symbols::SymbolicWorker& getSymbolicWorker() const;
    const symbols::NumericalWorker& getNumericalWorker() const;
    symbols::NumericalWorker& getNumericalWorker();

    bool is_s_squared_initialized() const;
    bool is_g_sz_squared_initialized() const;
    bool is_g_sz_squared_derivatives_initialized() const;
    bool is_isotropic_exchange_derivatives_initialized() const;
    bool is_zero_field_splitting_initialized() const;

    std::optional<std::reference_wrapper<const operators::Operator>>
        getOperator(common::QuantityEnum) const;
    std::optional<std::reference_wrapper<const operators::Operator>>
    getOperatorDerivative(common::QuantityEnum, const symbols::SymbolName&) const;
    const std::map<std::pair<common::QuantityEnum, symbols::SymbolName>, operators::Operator>&
    getOperatorDerivatives() const;
    const lexicographic::IndexConverter& getIndexConverter() const;

  private:
    symbols::NumericalWorker numericalWorker_;
    const lexicographic::IndexConverter converter_;
    struct OperatorsHistory {
        bool isotropic_exchange_in_hamiltonian = false;
        bool s_squared = false;
        bool g_sz_squared = false;
        bool isotropic_exchange_derivatives = false;
        bool g_sz_squared_derivatives = false;
        bool zfs_in_hamiltonian = false;
        bool zfs_derivative = false;
    };
    OperatorsHistory operators_history_;

    void InitializeIsotropicExchange();
    void InitializeIsotropicExchangeDerivatives();
    void InitializeZeroFieldSplitting();
    void InitializeZeroFieldSplittingDerivative();
    void InitializeSSquared();
    void InitializeGSzSquared();
    void InitializeGSzSquaredDerivatives();

    // std::shared_ptr<BasicQuantity>, where BasicQuantity is virtual class?
    // or Quantity can consist std::unique_ptr<VirtualMatrix>?
    // or can we just separate Operator and Spectrum?
    operators::Operator energy_operator;
    std::optional<operators::Operator> s_squared_operator;
    std::optional<operators::Operator> g_sz_squared_operator;
    std::map<std::pair<common::QuantityEnum, symbols::SymbolName>, operators::Operator>
        derivatives_map_;
};

}  // namespace model

#endif  //SPINNER_MODEL_H
