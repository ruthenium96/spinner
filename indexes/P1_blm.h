#ifndef JULY_P1_BLM_H
#define JULY_P1_BLM_H

#include "symm_blm.h"

// Character table for point group P1:
// P1 |  E
// -------
// A  | +1

struct Quantum_Numbers_P1 : Quantum_Numbers {
public:
    std::vector<unsigned int> path;

    Quantum_Numbers_P1(int total_mult_, int symmetry_, std::vector<unsigned int> path);

    // NB: decreasing order of total_mult!
    bool operator<(const Quantum_Numbers_P1 &other) const {
        return (this->total_mult > other.total_mult) ||
               ((this->total_mult == other.total_mult) && (this->path < other.path));
    }
};


struct Block_Lex_Map_P1 : Block_Lex_Map_Symm {
public:
    explicit Block_Lex_Map_P1(std::vector<int> mults);
    void construct_branching_diagram() override;

protected:
    double total_coeff(unsigned long state, unsigned long block) const override;

private:
    std::vector<Quantum_Numbers_P1> bds;
};

#endif