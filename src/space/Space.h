#ifndef JULY_SPACE_H
#define JULY_SPACE_H

#include <deque>
#include <iostream>
#include <map>

#include "Subspace.h"
#include "src/common/lexicographic/IndexConverter.h"

namespace space {
struct Space {
    explicit Space(uint32_t total_space_size);
    explicit Space(std::vector<Subspace>&& m);

    std::vector<Subspace> blocks;
};
}  // namespace space

std::ostream& operator<<(std::ostream& os, const space::Space& space);

#endif  // JULY_SPACE_H
