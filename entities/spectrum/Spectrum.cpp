#include "Spectrum.h"

std::ostream &operator<<(std::ostream &os, const Spectrum &spectrum) {
    for (const Subspectrum& subspectrum : spectrum.blocks) {
        os << subspectrum;
    }
    os << "------" << std::endl;
    return os;
}

Spectrum::Spectrum(std::vector<Subspectrum> &&m) : blocks(std::move(m)) {
}
