#include "Quantity.h"

std::string common::get_quantity_name(const QuantityEnum quantityEnum) {
    if (quantityEnum == QuantityEnum::Energy) {
        return "Energy";
    } else if (quantityEnum == QuantityEnum::S_total_squared) {
        return "S_total_squared";
    } else if (quantityEnum == QuantityEnum::gSz_total_squared) {
        return "gSz_total_squared";
    }
}