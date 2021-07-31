#ifndef JULY_BASE_BLM_C2_H
#define JULY_BASE_BLM_C2_H

#include <utility>
#include "base_blm.h"

// Character table for point group P2:
// P2 | (12) | (21) |
// ---|------|------|
// A  |  +1  |  +1  |
// B  |  +1  |  -1  |

// There is two possible types of orbits:
// (ab) : (12)(ab) = (ab); (21)(ab) = (ba)
// (aa) : (12)(aa) = (aa); (21)(aa) = (aa)


// If we have n mults and k (2k <= n) pairs of symmetry-bounded mults,
// than at first, we add mults in the pairs
// 1-1) (S_1 + S_2)
// 1-2) (S_3 + S_4)
// 1-i) (S_{2i-1} + S{2i})
// at second, step-by-step add the results of the previous step
// 2-1) (S_1 + S_2) + (S_3 + S_4)
// 2-2) ((S_1 + S_2) + (S_3 + S_4)) + (S_5 + S_6)
// 2-i) ((S_1 + S_2) + ... + (S_{2i-1} + S_{2i})) + (S_{2i+1} + S_{2i+2})
// and at the end we add (n - 2k) non-paired mults
// 3-1) (S_1 + ... S_{2k}) + S_{2k+1}
// 3-i) (S_1 + ... S_{2k+i-1}) + S_{2k+i}
class Spin_P2_Diagramm {
public:
    int total_mult;
    unsigned int representation;
    // i-th number of the first k numbers is the total mults of i-th pair.
    std::vector<int> pairs_mults;
    // i-th number is changing of cumulative spin during i-th addition of total spin of i-th pair.
    // (i+k)-th number is changing of cumulative spin during i-th addition of non-paired i-th spin.
    std::vector<unsigned int> path;

    // NB: decreasing order of back_mult!
    bool operator<(const Spin_P2_Diagramm &other) const {
        return (this->representation < other.representation) ||
               ((this->representation == other.representation) && (this->total_mult > other.total_mult));
    }
};


struct Indexes_P2 : Indexes {
public:
    // num_of_repr = 2;
    // group_size = 2;
    Indexes_P2(std::vector<int> mults_, unsigned int pairs_);

    std::vector<unsigned int> dim_of_repr;
    void calculate_dim_of_repr();

    void construct_symm_decomposition();

    double total_coeff(unsigned long state, unsigned long sym) const;

    unsigned long symmetrized_block(unsigned long block);

    void construct_branching_diagram() override;

    std::vector<Spin_P2_Diagramm> bds;


private:
    unsigned int pairs;
    static unsigned int character_multiplication(unsigned int a, unsigned int b);
};



#endif
