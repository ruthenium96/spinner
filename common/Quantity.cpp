#include "Quantity.h"

std::string common::get_quantity_name(const QuantityEnum quantityEnum) {
    if (quantityEnum == QuantityEnum::Energy) {
        return "Energy";
    } else if (quantityEnum == QuantityEnum::S_total_squared) {
        return "S_total_squared";
    }
}

std::ostream& operator<<(std::ostream& os, const common::QuantityEnum& quantity_enum) {
    os << "Quantity type: " << common::get_quantity_name(quantity_enum) << std::endl;
    return os;
}