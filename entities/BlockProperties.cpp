#include "BlockProperties.h"

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

std::ostream &operator<<(std::ostream &os, const BlockProperties &properties) {
    os << "Total n-projection: ";
    if (properties.n_proj.has_value()) {
        os << properties.n_proj.value();
    } else {
        os << "none";
    }
    os << '\n'
    << "    dimensionality: " << properties.dimensionality << '\n'
    << "        degeneracy: " << properties.degeneracy << '\n'
    << "    representation: " << properties.get_representation_name()
    << std::endl;
    return os;
}
