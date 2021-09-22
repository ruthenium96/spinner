#include "Logger.h"

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

std::ostream &operator<<(std::ostream &os, const Subspace &subspace) {
    os << subspace.properties;
    for (auto& m: subspace.basis) {
        for (auto& d: m) {
            os << d.second << "*[" << d.first << "] ";
        }
        os << std::endl;
    }
    os << std::endl;
    return os;
}

std::ostream &operator<<(std::ostream &os, const Space &space) {
    for (const Subspace& subspace : space.blocks) {
        os << subspace;
    }
    os << "------" << std::endl;
    return os;
}
