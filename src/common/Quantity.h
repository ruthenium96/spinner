#ifndef SPINNER_QUANTITYENUM_H
#define SPINNER_QUANTITYENUM_H

#include <ostream>
#include <string>

#include "src/entities/matrix/Matrix.h"
#include "src/entities/spectrum/Spectrum.h"

namespace common {
enum QuantityEnum {
    Energy,
    S_total_squared,
};

std::string get_quantity_name(common::QuantityEnum);

struct Quantity {
    Matrix matrix_;
    Spectrum spectrum_;
};
}  // namespace common

std::ostream& operator<<(std::ostream& os, const common::QuantityEnum& quantity);

#endif  //SPINNER_QUANTITYENUM_H
