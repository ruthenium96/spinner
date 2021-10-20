#include "Submatrix.h"

std::ostream &operator<<(std::ostream &os, const Submatrix &submatrix) {
    os << submatrix.properties;
    os << submatrix.submatrix;
    os << std::endl;
    return os;
}
