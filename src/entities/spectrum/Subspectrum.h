#ifndef SPINNER_SUBSPECTRUM_H
#define SPINNER_SUBSPECTRUM_H

#include "src/entities/BlockProperties.h"
#include "src/linalg_structures/AbstractDenseVector.h"
#include "src/entities/matrix/Submatrix.h"

namespace spinner {

struct Subspectrum {
    BlockProperties properties;
    std::unique_ptr<linalg_structures::AbstractDenseVector> raw_data;

    Subspectrum() = default;

    Subspectrum(
        std::unique_ptr<linalg_structures::AbstractDenseVector> raw_data_,
        BlockProperties properties_);
};

} // namespace spinner

#endif  //SPINNER_SUBSPECTRUM_H
