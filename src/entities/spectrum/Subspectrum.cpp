#include "Subspectrum.h"

#include <utility>

std::ostream& operator<<(std::ostream& os, const Subspectrum& subspectrum) {
    os << subspectrum.properties << std::endl;
    subspectrum.raw_data->print(os);
    os << std::endl;
    return os;
}

Subspectrum::Subspectrum(
    std::unique_ptr<quantum::linear_algebra::AbstractDenseVector> raw_data_,
    BlockProperties properties_) {
    raw_data = std::move(raw_data_);
    properties = std::move(properties_);
}
