#ifndef SPINNER_SUBSPACE_H
#define SPINNER_SUBSPACE_H

#include <cstdint>
#include "src/entities/BlockProperties.h"
#include "src/linalg_structures/AbstractSparseSemiunitaryMatrix.h"

namespace spinner::space {
struct Subspace {
    BlockProperties properties;
    std::unique_ptr<linalg_structures::AbstractSparseSemiunitaryMatrix> decomposition;

    explicit Subspace(
        std::unique_ptr<linalg_structures::AbstractSparseSemiunitaryMatrix>&&);
            
    uint32_t size() const;
};
} // namespace spinner::space

#endif  // SPINNER_SUBSPACE_H
