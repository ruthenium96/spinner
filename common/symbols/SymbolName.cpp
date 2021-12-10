#include "SymbolName.h"

namespace symbols {
bool SymbolName::operator<(const SymbolName& rhs) const {
    return name_ < rhs.name_;
}
bool SymbolName::operator>(const SymbolName& rhs) const {
    return rhs < *this;
}
bool SymbolName::operator<=(const SymbolName& rhs) const {
    return !(rhs < *this);
}
bool SymbolName::operator>=(const SymbolName& rhs) const {
    return !(*this < rhs);
}
bool SymbolName::operator==(const SymbolName& rhs) const {
    return name_ == rhs.name_;
}
bool SymbolName::operator!=(const SymbolName& rhs) const {
    return !(rhs == *this);
}
const std::string& SymbolName::get_name() const {
    return name_;
}
SymbolName::SymbolName(std::string name) : name_(std::move(name)) {}
}  // namespace symbols