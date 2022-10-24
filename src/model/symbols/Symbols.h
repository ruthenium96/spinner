#ifndef SPINNER_SYMBOLS_H
#define SPINNER_SYMBOLS_H

#include <map>
#include <memory>
#include <optional>
#include <stdexcept>

#include "SymbolName.h"
#include "src/group/Group.h"
#include "src/model/NumericalParameters.h"

namespace model::symbols {

enum SymbolTypeEnum { not_specified, J, g_factor, Theta };

// TODO: Can we separate "model" part and "optimize-values" part?
class Symbols {
  public:
    explicit Symbols(size_t number_of_spins);

    bool isAllGFactorsEqual() const;
    bool isGFactorInitialized() const;
    bool isIsotropicExchangeInitialized() const;
    bool isThetaInitialized() const;

    std::shared_ptr<const TwoDNumericalParameters<double>>
    constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name);
    std::shared_ptr<const TwoDNumericalParameters<double>> getIsotropicExchangeParameters() const;
    std::shared_ptr<const OneDNumericalParameters<double>> getGFactorParameters() const;
    std::shared_ptr<const TwoDNumericalParameters<double>> getGGFactorProductParameters() const;
    std::shared_ptr<const double> getThetaParameter() const;

    SymbolName getIsotropicExchangeSymbolName(size_t i, size_t j) const;
    SymbolName getGFactorSymbolName(size_t i) const;
    SymbolName getThetaSymbolName() const;

    Symbols& assignSymbolToIsotropicExchange(
        const SymbolName& symbol_name,
        size_t center_a,
        size_t center_b);
    Symbols& assignSymbolToGFactor(const SymbolName& symbol_name, size_t center_a);
    Symbols& assignSymbolToTheta(const SymbolName& symbol_name);

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

    // TODO: Can we move all of this to SymbolName and rename it to Symbol?
    //  maybe except of value.
    struct SymbolData {
        double value;
        bool is_changeable;
        SymbolTypeEnum type_enum;
    };

    std::vector<std::vector<SymbolName>> symbolic_isotropic_exchanges_;
    std::vector<SymbolName> symbolic_g_factors_;
    std::optional<SymbolName> symbolic_Theta_;

    void updateIsotropicExchangeParameters();
    void updateGFactorParameters();
    void updateThetaParameter();

    std::shared_ptr<TwoDNumericalParameters<double>> numeric_isotropic_exchanges_;
    std::shared_ptr<OneDNumericalParameters<double>> numeric_g_factors_;
    std::shared_ptr<TwoDNumericalParameters<double>> numeric_g_g_factors_product_;
    std::shared_ptr<double> numeric_Theta_;

    std::vector<std::shared_ptr<TwoDNumericalParameters<double>>>
        numeric_isotropic_exchange_derivatives_;

    std::map<SymbolName, SymbolData> symbolsMap;
};

}  // namespace model::symbols
#endif  //SPINNER_SYMBOLS_H
