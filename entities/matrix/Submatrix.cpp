#include "Submatrix.h"

std::ostream &operator<<(std::ostream &os, const Submatrix &submatrix) {
    os << submatrix.properties;
    os << submatrix.raw_data;
    os << std::endl;
    return os;
}
