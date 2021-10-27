#include "Subspectrum.h"

std::ostream &operator<<(std::ostream &os, const Subspectrum &subspectrum) {
    os << "ENERGY" << std::endl  << subspectrum.energy_raw_data << std::endl;
    for (size_t i = 0; i < subspectrum.non_energy_raw_data.size(); ++i) {
        os << "VALUE " << i << std::endl << subspectrum.non_energy_raw_data[i] << std::endl;
    }
    return os;
}
