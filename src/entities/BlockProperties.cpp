#include "BlockProperties.h"

namespace spinner {

std::string BlockProperties::get_representation_name() const {
    if (representation.empty()) {
        return "none";
    }
    std::string s;
    for (auto repr : representation) {
        s += (std::to_string(repr));
    }
    return s;
}

} // namespace spinner