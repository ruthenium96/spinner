#ifndef JULY_BASE_BLM_C2_H
#define JULY_BASE_BLM_C2_H

//TODO: сперва реализовать C2, потом уже сделать C1.

#include <utility>
#include "base_blm.h"

// Character table for point group P2:
// C2 |  E | C2
// ------------
// A  | +1 | +1
// B  | +1 | -1


struct Indexes_P2 : Indexes {
public:
    // TODO: implement block -> [sym] and sym -> [block]

    Indexes_P2(std::vector<int> mults_, unsigned int pairs_);

    unsigned int num_of_repr = 2;
    unsigned int group_size = 2;

    std::vector<unsigned int> dim_of_repr;
    void calculate_dim_of_repr();
    void contruct_symm_decomposition();


    std::vector<std::vector<unsigned long>> sym_sum_boundaries;

    unsigned long symmetrized_block(unsigned long block);


private:
    unsigned int pairs;
    std::vector<unsigned int> counters;
};



#endif
