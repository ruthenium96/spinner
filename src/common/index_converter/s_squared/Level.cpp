#include "Level.h"
#include <stdexcept>

namespace index_converter::s_squared {

Level::Level(std::shared_ptr<std::vector<spin_algebra::Multiplicity>> initialMultiplicities,
        size_t number_of_summations) {
    initialMultiplicities_ = std::move(initialMultiplicities);
    intermediateMultiplicities_.resize(number_of_summations, 0);
}

std::shared_ptr<std::vector<spin_algebra::Multiplicity>> Level::getInitialMultiplicities() const {
    return initialMultiplicities_;
}

spin_algebra::Multiplicity Level::getMultiplicity(size_t number) const {
    size_t initial_multiplicities_size = initialMultiplicities_->size();
    size_t intermediate_multiplicities_size = intermediateMultiplicities_.size();
    if (number < initial_multiplicities_size) {
        return initialMultiplicities_->at(number);
    } else {
        return intermediateMultiplicities_.at(number - initial_multiplicities_size);
    }
}

spin_algebra::Multiplicity Level::total() const {
    return intermediateMultiplicities_.back();
}

size_t Level::getSize() const {
    return initialMultiplicities_->size() + intermediateMultiplicities_.size();
}

void Level::setMultiplicity(size_t number, spin_algebra::Multiplicity multiplicity) {
    size_t initial_multiplicities_size = initialMultiplicities_->size();
    if (number < initial_multiplicities_size) {
        throw std::invalid_argument("Cannot set initial multiplicity");
    } else {
        size_t position_to_change = number - initial_multiplicities_size;
        if (intermediateMultiplicities_.at(position_to_change) == 0) {
            intermediateMultiplicities_[position_to_change] = multiplicity;
        } else {
            throw std::invalid_argument("Position was already set");
        }
    }
}

std::strong_ordering Level::operator<=>(const Level& rhs) const {
    for (int i = intermediateMultiplicities_.size() - 1; i >= 0; --i) {
        if (intermediateMultiplicities_[i] < rhs.intermediateMultiplicities_[i]) {
            return std::strong_ordering::less;
        }
        if (intermediateMultiplicities_[i] > rhs.intermediateMultiplicities_[i]) {
            return std::strong_ordering::greater;
        }
    }
    return std::strong_ordering::equal;
}

double Level::getSpin(size_t number) const {
    return ((double)getMultiplicity(number) - 1.0) / 2.0;
}

} // namespace index_converter::s_squared