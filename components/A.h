#ifndef JULY_A_H
#define JULY_A_H

#include <deque>
#include "B.h"
#include <iostream>

struct Task {
    std::vector<int> mults;

    bool is_Tz_sorted = false;
    bool is_C2_symmetrized = false;

    std::deque<Subspace> blocks;

    explicit Task(std::vector<int> mults_) : mults(mults_) {

        int tensor_size = 1;
        for (int mult : mults) {
            tensor_size *= mult;
        }

        blocks.emplace_back();
        Subspace& lex_block = blocks[0];
        lex_block.basis.resize(tensor_size);
        for (unsigned long i = 0; i < tensor_size; ++i) {
            lex_block.basis[i].push_back({i, 1.0});
        }
    }

    void print() {
        for (Subspace& Ss : blocks) {
            for (auto& v: Ss.basis) {
                for (Decomposition_& d: v) {
                    std::cout << d.coeff << "*" << d.index << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

};


#endif //JULY_A_H
