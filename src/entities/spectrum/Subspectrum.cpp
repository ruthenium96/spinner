#include "Subspectrum.h"

#include <utility>

namespace spinner {

Subspectrum::Subspectrum(
    std::unique_ptr<linalg_structures::AbstractDenseVector> raw_data_,
    BlockProperties properties_) {
    raw_data = std::move(raw_data_);
    properties = std::move(properties_);
}

} // namespace spinner