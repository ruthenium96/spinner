#include "ModelInput.h"

#include <cassert>

namespace model {

ModelInput::ModelInput(const std::vector<spin_algebra::Multiplicity>& mults) :
    symbols_(mults.size()),
    mults_(mults) {}

const symbols::SymbolicWorker& ModelInput::getSymbolicWorker() const {
    return symbols_;
}

symbols::SymbolicWorker& ModelInput::modifySymbolicWorker() {
    return symbols_;
}
const std::vector<spin_algebra::Multiplicity>& ModelInput::getMults() const {
    return mults_;
}

}  // namespace model