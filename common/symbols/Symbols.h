#ifndef JULY_SYMBOLS_H
#define JULY_SYMBOLS_H

#include <map>
#include <memory>
#include <stdexcept>

#include "SymbolName.h"
#include "entities/data_structures/DenseMatrix.h"
#include "group/Group.h"

namespace symbols {

enum SymbolTypeEnum { not_specified, J, g_factor };

class Symbols {
  public:
    explicit Symbols(size_t number_of_spins);

    [[nodiscard]] bool isAllGFactorsEqual() const;
    [[nodiscard]] bool symmetry_consistence(const group::Group& group) const;

    [[nodiscard]] std::shared_ptr<const DenseMatrix>
    constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name);
    [[nodiscard]] std::shared_ptr<const DenseMatrix> getIsotropicExchangeParameters() const;
    [[nodiscard]] std::shared_ptr<const DenseVector> getGFactorParameters() const;

    void assignSymbolToIsotropicExchange(
        const SymbolName& symbol_name,
        size_t center_a,
        size_t center_b);
    void assignSymbolToGFactor(const SymbolName& symbol_name, size_t center_a);

    SymbolName addSymbol(
        const std::string& name_string,
        double initial_value,
        bool is_changeable,
        SymbolTypeEnum type_enum);
    SymbolName addSymbol(const std::string& name, double initial_value, bool is_changeable);
    SymbolName addSymbol(const std::string& name, double initial_value);

    void setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value);

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

    std::vector<std::vector<SymbolName>> symbolic_isotropic_exchanges_;
    std::vector<SymbolName> symbolic_g_factors_;

    void updateIsotropicExchangeParameters();
    void updateGFactorParameters();

    std::shared_ptr<DenseMatrix> numeric_isotropic_exchanges_;
    std::shared_ptr<DenseVector> numeric_g_factors_;
    std::vector<std::shared_ptr<DenseMatrix>> numeric_isotropic_exchange_derivatives_;

    std::map<SymbolName, SymbolData> symbolsMap;
};

}  // namespace symbols
#endif  //JULY_SYMBOLS_H
