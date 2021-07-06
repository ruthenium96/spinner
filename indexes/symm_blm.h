#ifndef CJULY_SYMM_BLM_H
#define CJULY_SYMM_BLM_H

#include <utility>
#include "base_blm.h"
#include "../libs/wignerSymbols/include/wignerSymbols.h"
#include "boost/multi_array.hpp"

struct Quantum_Numbers {
public:
    int total_mult;
    int symmetry;
    Quantum_Numbers(int total_mult_, int symmetry_);
};

struct Spin_Boundary {
    int total_mult;
    std::vector<unsigned long> symm_boundaries;
};

struct Block_Lex_Map_Symm : Block_Lex_Map {
public:
    explicit Block_Lex_Map_Symm(std::vector<int> mults) : Block_Lex_Map(std::move(mults)) {}
    arma::dmat construct_transformation(int spin, int symmetry) const;
    std::vector<Spin_Boundary> spin_boundaries;

    template<typename QN_T>
    void counstruct_spin_boundaries(const std::vector<QN_T> &bds) {
        int current_spin = INT32_MAX;
        int current_symm = 0;

        for (unsigned long i = 0; i < bds.size(); ++i) {
            if (current_spin > bds[i].total_mult) {
                if (!spin_boundaries.empty()) {
                    spin_boundaries.back().symm_boundaries.push_back(i);
                }
                current_spin = bds[i].total_mult;
                current_symm = 0;
                spin_boundaries.push_back({current_spin, std::vector<unsigned long>(1, i)});
            }
            if (current_symm < bds[i].symmetry) {
                spin_boundaries.back().symm_boundaries.push_back(i);
                current_symm = bds[i].symmetry;
            }
        }
        spin_boundaries.back().symm_boundaries.push_back(bds.size());

//        for (const auto& v : spin_boundaries) {
//            for (auto n : v.symm_boundaries) {
//                std::cout << n << " ";
//            }
//            std::cout << std::endl;
//        }
    }

    // symmetry-dependent function:
    virtual void construct_branching_diagram() = 0;

protected:
    mutable boost::multi_array<double, 5> CGs;
    double max_t_spin;
    double max_n_spin;

    void resize_CGs();
    double hashed_clebsh_gordan(double l1, double l2, double l3, double m1, double m2) const;
    // symmetry-dependent functions:
    virtual double total_coeff(unsigned long state, unsigned long block) const = 0;

};

#endif
