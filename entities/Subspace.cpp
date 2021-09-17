#include "Subspace.h"

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
