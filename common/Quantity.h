#ifndef JULY_QUANTITYENUM_H
#define JULY_QUANTITYENUM_H

#include <string>
#include <ostream>

enum QuantityEnum {
    Energy,
    S_total_squared,
    };

[[nodiscard]] std::string get_quantity_name(QuantityEnum);


std::ostream &operator<<(std::ostream &os, const QuantityEnum &quantity);

#endif //JULY_QUANTITYENUM_H
