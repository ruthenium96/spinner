#ifndef SPINNER_SPECTRUM_H
#define SPINNER_SPECTRUM_H

#include <vector>

#include "Subspectrum.h"
#include "src/entities/matrix/Matrix.h"

struct Spectrum {
    Spectrum() = default;

    std::vector<Subspectrum> blocks;
    explicit Spectrum(std::vector<Subspectrum>&& m);
};

std::ostream& operator<<(std::ostream& os, const Spectrum& space);

#endif  //SPINNER_SPECTRUM_H
