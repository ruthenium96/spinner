#include "Quantity.h"

std::string get_quantity_name(const QuantityEnum quantityEnum) {
    if (quantityEnum == QuantityEnum::Energy) {
        return "Energy";
    } else if (quantityEnum == QuantityEnum::S_total_squared) {
        return "S_total_squared";
    }
}

std::ostream &operator<<(std::ostream &os, const QuantityEnum &quantity_enum) {
    os << "Quantity type: " << get_quantity_name(quantity_enum) << std::endl;
    return os;
}