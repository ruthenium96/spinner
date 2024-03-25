#ifndef SPINNER_MODELINPUT_H
#define SPINNER_MODELINPUT_H

#include <vector>

#include "src/common/Quantity.h"
#include "src/model/operators/Operator.h"
#include "src/model/symbols/SymbolicWorker.h"
namespace model {
// class ModelInput is responsible for user input: SymbolicWorker and Multiplicities
class ModelInput {
  public:
    explicit ModelInput(const std::vector<spin_algebra::Multiplicity>& mults);

    const symbols::SymbolicWorker& getSymbolicWorker() const;
    symbols::SymbolicWorker& modifySymbolicWorker();
    const std::vector<spin_algebra::Multiplicity>& getMults() const;
    bool operator==(const ModelInput& rhs) const = default;
    bool operator!=(const ModelInput& rhs) const = default;

  private:
    symbols::SymbolicWorker symbols_;
    std::vector<spin_algebra::Multiplicity> mults_;
};
}  // namespace model

#endif  //SPINNER_MODELINPUT_H
