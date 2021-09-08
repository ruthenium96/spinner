#ifndef JULY_SPACE_H
#define JULY_SPACE_H

#include "Entity.h"
#include "Subspace.h"
#include "common/LexicographicIndexConverter.h"
#include <deque>
#include <iostream>
#include <map>

struct Space : public entities::Entity {

    explicit Space(uint32_t total_space_size);
    friend std::ostream& operator<<(std::ostream& os, const Space& space);

    std::deque<Subspace> blocks;
};

#endif // JULY_SPACE_H
