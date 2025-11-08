#ifndef SPINNER_SUBSPACE_H
#define SPINNER_SUBSPACE_H

#include <cstdint>
#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"

namespace space {
struct Subspace {
    BlockProperties properties;
    std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix> decomposition;

    explicit Subspace(
        std::unique_ptr<quantum::linear_algebra::AbstractSparseSemiunitaryMatrix>&&);
            
    uint32_t size() const;
};
}  // namespace space

#endif  // SPINNER_SUBSPACE_H
