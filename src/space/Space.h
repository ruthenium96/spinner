#ifndef JULY_SPACE_H
#define JULY_SPACE_H

#include <deque>
#include <iostream>
#include <map>

#include "Subspace.h"
#include "src/common/lexicographic/IndexConverter.h"

namespace space {
// Space is responsible for _sparse_ unitary transformation of Hamiltonian matrix.
// This transformation leads to block-diagonal Hamiltonian matrix.
class Space {
  public:
    Space() = default;
    explicit Space(uint32_t total_space_size);
    explicit Space(std::vector<Subspace>&& m);

    std::vector<Subspace>& getBlocks();
    const std::vector<Subspace>& getBlocks() const;

  private:
    std::vector<Subspace> blocks_;
};
}  // namespace space

std::ostream& operator<<(std::ostream& os, const space::Space& space);

#endif  // JULY_SPACE_H
