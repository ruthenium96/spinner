#ifndef JULY_BASE_BLM_C1_H
#define JULY_BASE_BLM_C1_H

#include "base_blm.h"

// Character table for point group P1:
// P1 | (1) |
// ----------
// A  |  +1 |

// There is one possible type of orbit:
// (a) : (1)(a) = (a)

struct Indexes_P1 : Indexes {
    // num_of_repr = 1;
    // group_size = 1;

public:
    explicit Indexes_P1(std::vector<int> mults_);

    std::vector<unsigned int> dim_of_repr;

    void construct_branching_diagram() override;

    void construct_symm_decomposition();

};

#endif
