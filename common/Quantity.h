#ifndef JULY_QUANTITYENUM_H
#define JULY_QUANTITYENUM_H

#include <ostream>
#include <string>

#include "entities/matrix/Matrix.h"
#include "entities/operator/Operator.h"
#include "entities/spectrum/Spectrum.h"

namespace common {
enum QuantityEnum {
    Energy,
    S_total_squared,
};

[[nodiscard]] std::string get_quantity_name(common::QuantityEnum);

struct Quantity {
    Operator operator_;
    Matrix matrix_;
    Spectrum spectrum_;
};
}  // namespace common

std::ostream& operator<<(std::ostream& os, const common::QuantityEnum& quantity);

#endif  //JULY_QUANTITYENUM_H
