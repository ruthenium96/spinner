#include "P1_blm.h"

#include <utility>

Block_Lex_Map_P1::Block_Lex_Map_P1(std::vector<int> mults) : Block_Lex_Map_Symm(std::move(mults)) {
    max_t_spin = std::accumulate(spins.begin(), spins.end(), 0.0);
    max_n_spin = *std::max_element(spins.begin(), spins.end());
    resize_CGs();
}

// It is a branching diagram, i.e. the diagram produced by step-by-step addition of spins.
// 1) S_1
// 2) (S_1) + S_2
// 3) (S_1 + S_2) + S_3
// i) (S_1 + S_2 + ... + S_{i-1}) + S_i

// At the i-th place we have value of changing of cumulative spin during i-th addition.
void Block_Lex_Map_P1::construct_branching_diagram() {
    // TODO: rewrite without pushbacks!
    std::vector<Quantum_Numbers_P1> bd_from;
    std::vector<Quantum_Numbers_P1> bd_to;
    bd_from.emplace_back(1, 0, std::vector<unsigned int>{});

    for (int mult : mults) {
        bd_to.clear();
        for (const Quantum_Numbers_P1 & bdi : bd_from) {
            int t_mult = bdi.total_mult;
            std::vector<unsigned int> t_full = bdi.path;
            for (int k = std::abs(t_mult - mult) + 1; k < t_mult + mult; ++++k) {
                bd_to.emplace_back(k, 0, t_full);
                bd_to.back().path.push_back((k - t_mult + mult) / 2);
            }
        }
        std::swap(bd_from, bd_to);
    }

    std::stable_sort(bd_from.begin(), bd_from.end());
    bds = std::move(bd_from);

    counstruct_spin_boundaries(bds);

//    int current_spin = INT32_MAX;
//
//    for (unsigned long i = 0; i < bds.size(); ++i) {
//        if (current_spin > bds[i].total_mult) {
//            if (!spin_boundaries.empty()) {
//                spin_boundaries.back().symm_boundaries.push_back(i);
//            }
//            current_spin = bds[i].total_mult;
//            spin_boundaries.push_back({current_spin, std::vector<unsigned long>(1, i)});
//        }
//    }
//    spin_boundaries.back().symm_boundaries.push_back(bds.size());
}


double Block_Lex_Map_P1::total_coeff(unsigned long state, unsigned long block) const {
    // const std::vector<int> &nz_path
    const std::vector<unsigned int> & sp_path = bds[state].path;
    double c = 1;
    double t_spin = 0;
    double t_proj = 0;
    for (int i = 0; i < sp_path.size(); ++i) {
        double n_spin = ((double) mults[i] - 1.0) / 2.0;
        double n_proj = (double) block_to_nzi(block, i) - n_spin;
        double delta_t_spin = (double) sp_path[i] - n_spin;
        c *= hashed_clebsh_gordan(t_spin, n_spin,t_spin + delta_t_spin,
                                  t_proj, n_proj);
//            c *= WignerSymbols::clebschGordan(t_spin, n_spin,t_spin + delta_t_spin,
//                                              t_proj, n_proj,t_proj + n_proj);
        t_spin = t_spin + delta_t_spin;
        t_proj = t_proj + n_proj;
    }
    return c;
}

Quantum_Numbers_P1::Quantum_Numbers_P1(int total_mult_, int symmetry_, std::vector<unsigned int> path_) : Quantum_Numbers(
        total_mult_, symmetry_), path(std::move(path_)) {}
