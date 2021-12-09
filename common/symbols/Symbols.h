#ifndef JULY_SYMBOLS_H
#define JULY_SYMBOLS_H

#include <map>
#include <memory>
#include <stdexcept>

#include "entities/data_structures/DenseMatrix.h"
#include "group/Group.h"

namespace symbols {

enum SymbolTypeEnum { not_specified, J, g_factor };

class Symbols {
  public:
    struct SymbolName {
        bool operator<(const SymbolName& rhs) const;
        bool operator>(const SymbolName& rhs) const;
        bool operator<=(const SymbolName& rhs) const;
        bool operator>=(const SymbolName& rhs) const;
        bool operator==(const SymbolName& rhs) const;
        bool operator!=(const SymbolName& rhs) const;
        std::string name;
    };
    explicit Symbols(size_t number_of_spins);

    [[nodiscard]] bool hasIsotropicExchangeParameters() const;
    [[nodiscard]] bool isAllGFactorsEqual() const;
    [[nodiscard]] bool symmetry_consistence(const group::Group& group) const;

    [[nodiscard]] std::shared_ptr<const DenseMatrix>
    constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name);
    [[nodiscard]] std::shared_ptr<const DenseMatrix> constructIsotropicExchangeParameters();
    [[nodiscard]] std::shared_ptr<const DenseVector> constructGFactorParameters();

    void addIsotropicExchange(const SymbolName& symbol_name, size_t center_a, size_t center_b);
    void addGFactor(const SymbolName& symbol_name, size_t center_a);

    const SymbolName addSymbol(
        const std::string& name_string,
        double initial_value,
        bool is_changeable,
        SymbolTypeEnum type_enum);
    void setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value);
    void updateIsotropicExchangeParameters();
    void updateGFactorParameters();
    // TODO: implement it
    //    void updateIsotropicExchangeDerivativeParameters(const std::string& symbol_name);

    [[nodiscard]] std::vector<SymbolName> getChangeableNames(SymbolTypeEnum type_enum) const;
    [[nodiscard]] std::vector<SymbolName> getChangeableNames() const;
    [[nodiscard]] double getValueOfName(const SymbolName& symbol_name) const;

  private:
    size_t number_of_spins_;

    struct SymbolData {
        double value;
        bool is_changeable;
        SymbolTypeEnum type_enum;
    };

    std::vector<std::vector<SymbolName>> isotropic_exchange_parameters_names_;
    std::vector<SymbolName> g_factor_names_;

    std::shared_ptr<DenseMatrix> isotropic_exchange_parameters_values_;
    std::shared_ptr<DenseVector> g_factor_values_;
    std::vector<std::shared_ptr<DenseMatrix>> isotropic_exchange_derivatives_values_;

    std::map<SymbolName, SymbolData> symbols_;
};

}  // namespace symbols
#endif  //JULY_SYMBOLS_H
