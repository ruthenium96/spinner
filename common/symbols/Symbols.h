#ifndef JULY_SYMBOLS_H
#define JULY_SYMBOLS_H

#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "entities/data_structures/DenseMatrix.h"
#include "group/Group.h"

namespace symbols {

class Symbols {
  public:
    explicit Symbols(size_t number_of_spins);

    [[nodiscard]] bool hasIsotropicExchangeParameters() const;
    [[nodiscard]] std::shared_ptr<const DenseMatrix> constructIsotropicExchangeParameters();
    [[nodiscard]] std::shared_ptr<const DenseVector> constructGFactorParameters();

    void addIsotropicExchange(const std::string& symbol_name, size_t center_a, size_t center_b);
    void addGFactor(const std::string& symbol_name, size_t center_a);

    void addSymbol(const std::string& name, double initial_value, bool is_changeable);

    [[nodiscard]] bool symmetry_consistence(const group::Group& group) const;

  private:
    size_t number_of_spins_;

    std::vector<std::vector<std::string>> isotropic_exchange_parameters_symbols_;
    std::shared_ptr<DenseMatrix> isotropic_exchange_parameters_values_;

    std::vector<std::string> g_factor_symbols_;
    std::shared_ptr<DenseVector> g_factor_values_;

    std::unordered_map<std::string, double> name_to_value_;
    std::vector<std::string> changeable_symbols_;
};

}  // namespace symbols
#endif  //JULY_SYMBOLS_H
