#ifndef SPINNER_SUBSPACE_H
#define SPINNER_SUBSPACE_H

#include <cstdint>
#include "src/entities/BlockProperties.h"
#include "src/entities/data_structures/AbstractSparseSemiunitaryMatrix.h"

namespace spinner::space {
struct Subspace {
    BlockProperties properties;
    std::unique_ptr<linear_algebra::AbstractSparseSemiunitaryMatrix> decomposition;

    explicit Subspace(
        std::unique_ptr<linear_algebra::AbstractSparseSemiunitaryMatrix>&&);
            
    uint32_t size() const;
};
} // namespace spinner::space

#endif  // SPINNER_SUBSPACE_H
