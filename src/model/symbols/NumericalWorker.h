#ifndef SPINNER_NUMERICALWORKER_H
#define SPINNER_NUMERICALWORKER_H

#include "SymbolicWorker.h"

namespace model::symbols {

class NumericalWorker {
  public:
    NumericalWorker(SymbolicWorker symbols, size_t number_of_spins);
    const SymbolicWorker& getSymbolicWorker() const;
    void setNewValueToChangeableSymbol(const SymbolName& symbol_name, double new_value);

    std::shared_ptr<const TwoDNumericalParameters<double>>
    constructIsotropicExchangeDerivativeParameters(const SymbolName& symbol_name);
    std::shared_ptr<const OneDNumericalParameters<double>>
    constructZFSDerivativeParameters(const SymbolName& symbol_name);
    std::pair<
        std::shared_ptr<const OneDNumericalParameters<double>>,
        std::shared_ptr<const TwoDNumericalParameters<double>>>
    constructGGDerivativeParameters(const SymbolName& symbol_name);

    std::shared_ptr<const TwoDNumericalParameters<double>> getIsotropicExchangeParameters() const;
    std::shared_ptr<const OneDNumericalParameters<double>> getGFactorParameters() const;
    std::pair<
        std::shared_ptr<const OneDNumericalParameters<double>>,
        std::shared_ptr<const TwoDNumericalParameters<double>>>
    getGGParameters() const;
    std::shared_ptr<const double> getThetaParameter() const;
    std::pair<
        std::shared_ptr<OneDNumericalParameters<double>>,
        std::optional<std::shared_ptr<OneDNumericalParameters<double>>>>
    getZFSParameters() const;

  private:
    void updateIsotropicExchangeParameters();
    void updateGFactorParameters();
    void updateZFSParameters();
    void updateGGFactorDerivativesParameters(const SymbolName& symbol_name);
    void updateThetaParameter();

    std::shared_ptr<TwoDNumericalParameters<double>> numeric_isotropic_exchanges_;
    std::shared_ptr<OneDNumericalParameters<double>> numeric_g_factors_;
    std::pair<
        std::shared_ptr<OneDNumericalParameters<double>>,
        std::optional<std::shared_ptr<OneDNumericalParameters<double>>>>
        numeric_ZFS_;
    std::shared_ptr<double> numeric_Theta_;

    std::pair<
        std::shared_ptr<OneDNumericalParameters<double>>,
        std::shared_ptr<TwoDNumericalParameters<double>>>
        numeric_g_g_;
    std::map<
        SymbolName,
        std::pair<
            std::shared_ptr<OneDNumericalParameters<double>>,
            std::shared_ptr<TwoDNumericalParameters<double>>>>
        numeric_g_g_derivatives_;

    std::vector<std::shared_ptr<TwoDNumericalParameters<double>>>
        numeric_isotropic_exchange_derivatives_;

    SymbolicWorker symbolicWorker_;
    size_t number_of_spins_;
};

}  // namespace model::symbols

#endif  //SPINNER_NUMERICALWORKER_H
