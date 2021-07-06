#include "P2_blm.h"

#include <utility>

Block_Lex_Map_P2::Block_Lex_Map_P2(std::vector<int> mults,
                                   unsigned int pairs_) : Block_Lex_Map_Symm(std::move(mults)),
                                   pairs(pairs_) {
    max_t_spin = std::accumulate(spins.begin(), spins.end(), 0.0);
    max_n_spin = 2 * (*std::max_element(spins.begin(), spins.end()));
    resize_CGs();
}

void Block_Lex_Map_P2::construct_branching_diagram() {
    std::vector<Spin_P2_Diagramm> bd_from;
    std::vector<Spin_P2_Diagramm> bd_to;
    bd_from.push_back({1, 0});

    // сложение внутри пар:
    for (int i = 0; i < pairs; ++i) {
        bd_to.clear();
        for (const Spin_P2_Diagramm & bdi : bd_from) {
            const std::vector<int>& pm_full = bdi.pairs_mults;
            unsigned int symm_of_new_pair = 0;
            for (int m = 2 * mults[2 * i] - 1; m > 0 ; ----m) {
                unsigned int cum_symm = character_multiplication(symm_of_new_pair, bdi.symmetry);
                bd_to.push_back({1, cum_symm, pm_full, std::vector<unsigned int>()});
                bd_to.back().pairs_mults.push_back(m);
                // TODO: effective, but dirty:
                symm_of_new_pair = 1 - symm_of_new_pair;
            }
        }
        std::swap(bd_from, bd_to);
    }

    const std::vector<Spin_P2_Diagramm> bd_copy = bd_from;

    std::vector<Spin_P2_Diagramm> bd_result;

    for (const Spin_P2_Diagramm & sp2d : bd_copy) {
        const std::vector<int> & v = sp2d.pairs_mults;
        bd_from.clear();
        bd_from.push_back({1, sp2d.symmetry, v, std::vector<unsigned int>{}});
        for (int mult : v) {
            bd_to.clear();
            for (const Spin_P2_Diagramm & bdi : bd_from) {
                int t_mult = bdi.total_mult;
                std::vector<unsigned int> t_full = bdi.path;
                for (int k = std::abs(t_mult - mult) + 1; k < t_mult + mult; ++++k) {
                    bd_to.push_back({k, bdi.symmetry, bdi.pairs_mults, t_full});
                    bd_to.back().path.push_back((k - t_mult + mult) / 2);
                }
            }
            std::swap(bd_from, bd_to);
        }
        for (unsigned int i = 2 * pairs; i < mults.size(); ++i) {
            int mult = mults[i];
            bd_to.clear();
            for (const Spin_P2_Diagramm & bdi : bd_from) {
                int t_mult = bdi.total_mult;
                std::vector<unsigned int> t_full = bdi.path;
                for (int k = std::abs(t_mult - mult) + 1; k < t_mult + mult; ++++k) {
                    bd_to.push_back({k, bdi.symmetry, bdi.pairs_mults, t_full});
                    bd_to.back().path.push_back((k - t_mult + mult) / 2);
                }
            }
            std::swap(bd_from, bd_to);
        }
        bd_result.insert(bd_result.end(), bd_from.begin(), bd_from.end());
    }

    std::stable_sort(bd_result.begin(), bd_result.end());
    bds = std::move(bd_result);

//    for (const Spin_P2_Diagramm & pmf : bds) {
//        for (auto n : pmf.pairs_mults) {
//            std::cout << n << " ";
//        }
//        std::cout << "| ";
//        for (auto n : pmf.path) {
//            std::cout << n << " ";
//        }
//        std::cout << "| symmetry : " << pmf.symmetry << " | total_mult: " << pmf.total_mult << std::endl;
//    }
//    std::cout << std::endl;

    counstruct_spin_boundaries(bds);
}

unsigned int Block_Lex_Map_P2::character_multiplication(unsigned int a, unsigned int b) {
    // A * A == A => (0 + 0) == 0 (mod 2)
    // A * B == B => (0 + 1) == 1 (mod 2)
    // B * B == A => (1 + 1) == 0 (mod 2)
    return (a + b) % 2;
}

double Block_Lex_Map_P2::total_coeff(unsigned long state, unsigned long block) const {
//    const std::vector<int> &nz_path
    const std::vector<int> & pr_mult = bds[state].pairs_mults;
    const std::vector<unsigned int> & sp_path = bds[state].path;
    double c = 1;
    double t_spin = 0;
    double t_proj = 0;
    for (unsigned int i = 0; i < pairs; ++i) {
        // pair part
        double f_spin = ((double) mults[2 * i] - 1.0) / 2.0; // first spin in tha pair
        double s_spin = f_spin;                              // second spin in the pair
        double p_spin = ((double) pr_mult[i] - 1.0) / 2.0;   // total spin of the pair
        double f_proj = (double) block_to_nzi(block, 2 * i) - f_spin;
        double s_proj = (double) block_to_nzi(block, 2 * i + 1) - s_spin;
        c *= hashed_clebsh_gordan(f_spin, s_spin, p_spin,
                                  f_proj, s_proj);

        // total part
        double p_proj = f_proj + s_proj;
        double delta_t_spin = (double) sp_path[i] - p_spin;
        c *= hashed_clebsh_gordan(t_spin, p_spin,t_spin + delta_t_spin,
                                  t_proj, p_proj);
        t_spin = t_spin + delta_t_spin;
        t_proj = t_proj + p_proj;
    }
    for (unsigned int j = 2 * pairs; j < sp_path.size(); ++j) {
        double n_spin = ((double) mults[j] - 1.0) / 2.0;
        double n_proj = (double) block_to_nzi(block, j) - n_spin;
        double delta_t_spin = (double) sp_path[j] - n_spin;
        c *= hashed_clebsh_gordan(t_spin, n_spin,t_spin + delta_t_spin,
                                  t_proj, n_proj);
        t_spin = t_spin + delta_t_spin;
        t_proj = t_proj + n_proj;
    }
    return c;
}

