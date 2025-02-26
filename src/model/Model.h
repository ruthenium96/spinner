#ifndef SPINNER_MODEL_H
#define SPINNER_MODEL_H

#include <memory>
#include "ModelInput.h"
#include "symbols/NumericalWorker.h"

#include "src/model/operators/AbstractOperatorConstructor.h"

namespace model {

class Model {
  public:
    explicit Model(
        ModelInput modelInput, 
        std::unique_ptr<operators::AbstractOperatorConstructor>&& operator_constructor);

    // TODO: remove it?
    const symbols::SymbolicWorker& getSymbolicWorker() const;
    const symbols::NumericalWorker& getNumericalWorker() const;

    void
    setNewValueToChangeableSymbol(const model::symbols::SymbolName& symbol_name, double new_value);

    bool is_g_sz_squared_initialized() const;
    bool is_g_sz_squared_derivatives_initialized() const;

    bool gFactorsAreAllNoneOrAreTheSame() const;

    std::optional<std::shared_ptr<const operators::Operator>>
        getOperator(common::QuantityEnum) const;
    std::optional<std::shared_ptr<const operators::Operator>>
    getOperatorDerivative(common::QuantityEnum, const symbols::SymbolName&) const;

    // These methods are not 'const' because they return non-const pointers.
    const std::map<common::QuantityEnum, std::shared_ptr<operators::Operator>>& getOperators();
    const std::map<
        std::pair<common::QuantityEnum, symbols::SymbolName>,
        std::shared_ptr<operators::Operator>>&
    getOperatorDerivatives();

  private:
    symbols::NumericalWorker numericalWorker_;
    std::unique_ptr<operators::AbstractOperatorConstructor> operator_constructor_;

    struct OperatorsHistory {
        bool isotropic_exchange_in_hamiltonian = false;
        bool s_squared = false;
        bool g_sz_squared = false;
        bool isotropic_exchange_derivatives = false;
        bool g_sz_squared_derivatives = false;
        bool zfs_in_hamiltonian = false;
        bool zfs_derivative = false;
        bool isotropic_exchange_in_ito_hamiltonian = false;
    };
    OperatorsHistory operators_history_;

    symbols::NumericalWorker& getNumericalWorker();

    std::map<common::QuantityEnum, std::shared_ptr<operators::Operator>> operators_map_;
    std::map<
        std::pair<common::QuantityEnum, symbols::SymbolName>,
        std::shared_ptr<operators::Operator>>
        derivatives_map_;
};

}  // namespace model

#endif  //SPINNER_MODEL_H
