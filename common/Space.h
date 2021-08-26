#ifndef JULY_SPACE_H
#define JULY_SPACE_H

#include <deque>
#include <map>
#include "Subspace.h"
#include <iostream>

struct Space {
    std::vector<int> mults;

    bool is_Tz_sorted = false;
    bool is_C2_symmetrized = false;

    std::deque<Subspace> blocks;

    explicit Space(std::vector<int> mults_) : mults(mults_) {

        int tensor_size = 1;
        for (int mult : mults) {
            tensor_size *= mult;
        }

        blocks.emplace_back();
        Subspace& lex_block = blocks[0];
        lex_block.basis.resize(tensor_size);
        for (unsigned long i = 0; i < tensor_size; ++i) {
            lex_block.basis[i][i] = 1.0;
        }
    }

    void print() {
        for (Subspace& Ss : blocks) {
            for (auto& m: Ss.basis) {
                for (auto& d: m) {
                    std::cout << d.second << "*[" << d.first << "] ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << "------" << std::endl;
    }

};


#endif //JULY_SPACE_H
