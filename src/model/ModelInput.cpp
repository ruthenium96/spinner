#include "ModelInput.h"

#include <cassert>

namespace model {

ModelInput::ModelInput(const std::vector<int>& mults) : symbols_(mults.size()), mults_(mults) {}

const symbols::Symbols& ModelInput::getSymbols() const {
    return symbols_;
}

symbols::Symbols& ModelInput::getSymbols() {
    return symbols_;
}
const std::vector<int>& ModelInput::getMults() const {
    return mults_;
}

}  // namespace model