#ifndef SPINNER_SPECTRUM_H
#define SPINNER_SPECTRUM_H

#include <functional>
#include <vector>

#include "Subspectrum.h"

struct Spectrum {
    Spectrum() = default;

    std::vector<Subspectrum> blocks;
    explicit Spectrum(std::vector<Subspectrum>&& m);
};

struct SpectrumRef {
    std::vector<std::reference_wrapper<const Subspectrum>> blocks;

    explicit SpectrumRef(const Spectrum& spectrum);
};

#endif  //SPINNER_SPECTRUM_H
