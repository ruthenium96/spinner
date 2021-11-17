#ifndef JULY_SYMBOLS_H
#define JULY_SYMBOLS_H

#include <memory>
#include <stdexcept>

#include "entities/data_structures/DenseMatrix.h"

namespace symbols {

class Symbols {
  public:
    explicit Symbols(size_t number_of_spins);

    [[nodiscard]] bool hasIsotropicExchangeParameters() const;
    [[nodiscard]] std::shared_ptr<const DenseMatrix> getIsotropicExchangeParameters() const;

    void addIsotropicExchange(double value, size_t center_a, size_t center_b);
    void addIsotropicExchangeMatrix(DenseMatrix);

    //    enum SymbolTypeEnum {
    //        g,
    //        J,
    //    };
    //
    //    struct Symbol {
    //        std::string name;
    //        SymbolTypeEnum type;
    //        bool is_changeable;
    //    };
  private:
    size_t number_of_spins_;
    std::shared_ptr<DenseMatrix> isotropic_exchange_parameters_;
};

}  // namespace symbols
#endif  //JULY_SYMBOLS_H
