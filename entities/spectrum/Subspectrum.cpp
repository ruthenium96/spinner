#include "Subspectrum.h"

std::ostream &operator<<(std::ostream &os, const Subspectrum &subspectrum) {
    os << subspectrum.properties << std::endl;
    os << subspectrum.raw_data << std::endl;
    return os;
}
