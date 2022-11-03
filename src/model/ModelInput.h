#ifndef SPINNER_MODELINPUT_H
#define SPINNER_MODELINPUT_H

#include <vector>

#include "src/common/Quantity.h"
#include "src/model/operators/Operator.h"
#include "src/model/symbols/Symbols.h"
namespace model {
// class ModelInput is responsible for user input: Symbols and Multiplicities
class ModelInput {
  public:
    explicit ModelInput(const std::vector<int>& mults);

    const symbols::Symbols& getSymbols() const;
    symbols::Symbols& getSymbols();
    const std::vector<int>& getMults() const;

  private:
    symbols::Symbols symbols_;
    std::vector<int> mults_;
};
}  // namespace model

#endif  //SPINNER_MODELINPUT_H
