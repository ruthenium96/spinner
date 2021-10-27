#ifndef JULY_SPECTRUM_H
#define JULY_SPECTRUM_H

#include "Subspectrum.h"

#include <vector>


struct Spectrum {

    Spectrum() = default;

    std::vector<Subspectrum> blocks;
    explicit Spectrum(std::vector<Subspectrum>&& m);
};

std::ostream &operator<<(std::ostream &os, const Spectrum &space);

#endif //JULY_SPECTRUM_H
