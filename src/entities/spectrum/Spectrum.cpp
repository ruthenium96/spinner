#include "Spectrum.h"

Spectrum::Spectrum(std::vector<Subspectrum>&& m) : blocks(std::move(m)) {}

SpectrumRef::SpectrumRef(const Spectrum& spectrum) {
    for (int i = 0; i < spectrum.blocks.size(); ++i) {
        blocks.push_back(spectrum.blocks[i]);
    }
}
