#ifndef JULY_BASE_BLM_C1_H
#define JULY_BASE_BLM_C1_H

#include "base_blm.h"

// Character table for point group P1:
// P1 | (1) |
// ----------
// A  |  +1 |

// There is one possible type of orbit:
// (a) : (1)(a) = (a)

// It is a "branching diagram", i.e. the diagram produced by step-by-step addition of mults.
// 1) S_1
// 2) (S_1) + S_2
// 3) (S_1 + S_2) + S_3
// i) (S_1 + S_2 + ... + S_{i-1}) + S_i

struct Spin_P1_Diagramm {
public:
    int total_mult;
    // At the i-th place we have value of changing of cumulative spin during i-th addition.
    std::vector<unsigned int> path;

    // NB: decreasing order of back_mult!
    bool operator<(const Spin_P1_Diagramm &other) const {
        return (this->total_mult > other.total_mult);
    }
};


struct Indexes_P1 : Indexes {
    // num_of_repr = 1;
    // group_size = 1;

public:
    explicit Indexes_P1(std::vector<int> mults_);

    std::vector<unsigned int> dim_of_repr;

    void construct_branching_diagram() override;

    void construct_symm_decomposition();


private:
    std::vector<Quantum_Numbers> qns;

};

#endif
