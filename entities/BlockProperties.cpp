#include "BlockProperties.h"

std::ostream &operator<<(std::ostream &os, const BlockProperties &properties) {
    os << "Total n-projection: " << properties.n_proj << '\n'
       << "    dimensionality: " << properties.dimensionality << '\n'
       << "        degeneracy: " << properties.degeneracy << '\n'
       << "    representation: " << properties.get_representation_name()
       << std::endl;
    return os;
}

std::string BlockProperties::get_representation_name() const {
    std::string s = "";
    for (auto repr : representation) {
        s += (repr + '0');
    }
    return std::move(s);
}
