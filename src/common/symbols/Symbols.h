#ifndef JULY_SYMBOLS_H
#define JULY_SYMBOLS_H

#include <map>
#include <memory>
#include <stdexcept>

#include "SymbolName.h"
#include "src/entities/data_structures/DenseMatrix.h"
#include "src/group/Group.h"

namespace symbols {

enum SymbolTypeEnum { not_specified, J, g_factor };

class Symbols {
  public:
    explicit Symbols(size_t number_of_spins);

    bool isAllGFactorsEqual() const;
    bool isGFactorInitialized() const;
    bool isIsotropicExchangeInitialized() const;
    bool symmetry_consistence(const group::Group& group) const;

    std::shared_ptr<const DenseMatrix>
    constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name);
    std::shared_ptr<const DenseMatrix> getIsotropicExchangeParameters() const;
    std::shared_ptr<const DenseVector> getGFactorParameters() const;

    Symbols& assignSymbolToIsotropicExchange(
        const SymbolName& symbol_name,
        size_t center_a,
        size_t center_b);
    Symbols& assignSymbolToGFactor(const SymbolName& symbol_name, size_t center_a);

    SymbolName addSymbol(
        const std::string& name_string,
        double initial_value,
        bool is_changeable,
        SymbolTypeEnum type_enum);
    SymbolName addSymbol(const std::string& name, double initial_value, bool is_changeable);
    SymbolName addSymbol(const std::string& name, double initial_value);

    Symbols& setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value);

    std::vector<SymbolName> getChangeableNames(SymbolTypeEnum type_enum) const;
    std::vector<SymbolName> getChangeableNames() const;
    double getValueOfName(const SymbolName& symbol_name) const;

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
