#ifndef JULY_P2_BLM_H
#define JULY_P2_BLM_H

#include "symm_blm.h"

// Character table for point group P2:
// P2 |  E | C2
// ------------
// A  | +1 | +1
// B  | +1 | -1


// If we have n spins and k (2k <= n) pairs of symmetry-bounded spins,
// than at first, we add spins in the pairs
// 1-1) (S_1 + S_2)
// 1-2) (S_3 + S_4)
// 1-i) (S_{2i-1} + S{2i})
// at second, step-by-step add the results of the previous step
// 2-1) (S_1 + S_2) + (S_3 + S_4)
// 2-2) ((S_1 + S_2) + (S_3 + S_4)) + (S_5 + S_6)
// 2-i) ((S_1 + S_2) + ... + (S_{2i-1} + S_{2i})) + (S_{2i+1} + S_{2i+2})
// and at the end we add (n - 2k) non-paired spins
// 3-1) (S_1 + ... S_{2k}) + S_{2k+1}
// 3-i) (S_1 + ... S_{2k+i-1}) + S_{2k+i}

class Spin_P2_Diagramm {
public:
    int total_mult;
    unsigned int symmetry;
    // i-th number of the first k numbers is the total spins of i-th pair.
    std::vector<int> pairs_mults;
    // i-th number is changing of cumulative spin during i-th addition of total spin of i-th pair.
    // (i+k)-th number is changing of cumulative spin during i-th addition of non-paired i-th spin.
    std::vector<unsigned int> path;

    // NB: decreasing order of total_mult!
    bool operator<(const Spin_P2_Diagramm &other) const {
        return (this->total_mult > other.total_mult) ||
               ((this->total_mult == other.total_mult) && (this->symmetry < other.symmetry));
    }
};


struct Block_Lex_Map_P2 : Block_Lex_Map_Symm {
public:
    explicit Block_Lex_Map_P2(std::vector<int> mults, unsigned int pairs_);
    void construct_branching_diagram() override;

private:
    const unsigned int pairs;
    std::vector<Spin_P2_Diagramm> bds;
    // TODO: сделать pure virtual в классе-
    // TODO: имплементация, увы, только для абелевых групп.
    unsigned int character_multiplication(unsigned int a, unsigned int b);
    double total_coeff(unsigned long state, unsigned long block) const override;

};

#endif
