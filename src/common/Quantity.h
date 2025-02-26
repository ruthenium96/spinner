#ifndef SPINNER_QUANTITYENUM_H
#define SPINNER_QUANTITYENUM_H

#include "src/entities/matrix/Matrix.h"
#include "src/entities/spectrum/Spectrum.h"

namespace common {
enum QuantityEnum {
    Energy,
    S_total_squared, M_total_squared,
    g_squared_T00, gSz_total_squared };

struct Quantity {
    Matrix matrix_;
    Spectrum spectrum_;
};
}  // namespace common

#endif  //SPINNER_QUANTITYENUM_H
