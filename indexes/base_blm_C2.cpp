#include "base_blm_C2.h"

Indexes_P2::Indexes_P2(std::vector<int> mults_,
                       unsigned int pairs_) : Indexes(std::move(mults_)), pairs(pairs_) {
    num_of_repr = 2;
    group_size = 2;

    dim_of_repr.resize(num_of_repr);
    sym_sum_boundaries.resize(num_of_repr);

    block2sym.resize(tensor_size);
    sym2block.resize(tensor_size);
    calculate_dim_of_repr();
    // TODO: абстрагировать строки ниже
    construct_symm_decomposition();

//    for (const auto& vd : block2sym) {
//        for (auto d : vd) {
//            std::cout << "(" << d.coeff << " * |" << d.index << ">) + ";
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    for (const auto& vd : sym2block) {
//        for (auto d : vd) {
//            std::cout << "(" << d.coeff << " * |" << d.index << ">) + ";
//        }
//        std::cout << std::endl;
//    }
//
//    std::cout << std::endl;
//
//    for (const auto& repr : sym_sum_boundaries) {
//        for (auto n : repr) {
//            std::cout << n << " ";
//        }
//        std::cout << std::endl;
//    }

    // TODO: refactor it.
    max_t_spin = std::accumulate(spins.begin(), spins.end(), 0.0);
    max_n_spin = 2 * (*std::max_element(spins.begin(), spins.end()));
    resize_CGs();

}

// TODO: абстрагировать этот код!
void Indexes_P2::construct_symm_decomposition() {
    std::vector<unsigned int> counters;
    counters.resize(num_of_repr);
    counters[0] = 0;
    counters[1] = dim_of_repr[0];

    double sqrt_of_two = 1 / sqrt(2.0);
    unsigned long current_sum = UINT64_MAX;
    for (unsigned long block = 0; block < tensor_size; ++block) {
        std::vector<int> nzs = block_to_nzs(block);
        unsigned long sum = std::accumulate(nzs.begin(), nzs.end(), 0);
        if (sum != current_sum) {
            sym_sum_boundaries[0].push_back(counters[0]);
            sym_sum_boundaries[1].push_back(counters[1]);
            current_sum = sum;
        }
        if (block2sym[block].empty()) {
            unsigned long symm_block = symmetrized_block(block);
            if (symm_block == block) {
                block2sym[block] = std::vector<Decomposition>(1);
                block2sym[block][0] = {counters[0], 1.0};

                sym2block[counters[0]] = std::vector<Decomposition>(1);
                sym2block[counters[0]][0] = {block, 1.0};

                ++counters[0];
            } else {
                block2sym[block] = std::vector<Decomposition>(2);
                block2sym[symm_block] = std::vector<Decomposition>(2);
                block2sym[block][0] = {counters[0], sqrt_of_two};
                block2sym[block][1] = {counters[1], sqrt_of_two};
                block2sym[symm_block][0] = {counters[0], sqrt_of_two};
                block2sym[symm_block][1] = {counters[1], -sqrt_of_two};

                sym2block[counters[0]] = std::vector<Decomposition>(2);
                sym2block[counters[1]] = std::vector<Decomposition>(2);
                sym2block[counters[0]][0] = {block, sqrt_of_two};
                sym2block[counters[0]][1] = {symm_block, sqrt_of_two};
                sym2block[counters[1]][0] = {block, sqrt_of_two};
                sym2block[counters[1]][1] = {symm_block, -sqrt_of_two};

                ++counters[0];
                ++counters[1];
            }
        }
    }
    sym_sum_boundaries[0].push_back(sym_sum_boundaries[1][0]);
    sym_sum_boundaries[1].push_back(tensor_size);
}

void Indexes_P2::calculate_dim_of_repr() {
    unsigned int partitive_pair_product  = 1;
    for (int i = 0; i < 2 * pairs; ++++i) {
        partitive_pair_product *= mults[i];
    }

    unsigned int num_self_symmetrized = tensor_size / partitive_pair_product;
    unsigned int num_non_self_symmetrized = tensor_size - num_self_symmetrized;

    dim_of_repr[0] = num_self_symmetrized + num_non_self_symmetrized / 2;
    dim_of_repr[1] = num_non_self_symmetrized / 2;
}

unsigned long Indexes_P2::symmetrized_block(unsigned long block) {
    std::vector<int> nzs = block_to_nzs(block);
    for (int k = 0; k < pairs; ++k) {
        std::swap(nzs[2 * k], nzs[2 * k + 1]);
    }
    return nzs_to_block(nzs);
}

double Indexes_P2::total_coeff(unsigned long state, unsigned long sym) const {
    // first decomposition element always exists and its coefficient is non-zero in this group:
    unsigned long block = sym2block[sym][0].index;
    const std::vector<int> & pr_mult = bds[state].pairs_mults;
    const std::vector<unsigned int> & sp_path = bds[state].path;
    double c = 1;
    double t_spin = 0;
    double t_proj = 0;
    for (unsigned int i = 0; i < pairs; ++i) {
        // pair part
        double f_spin = ((double) mults[2 * i] - 1.0) / 2.0; // first spin in the pair
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

        if (c == 0.0) {
            return c;
        }
    }
    for (unsigned int j = pairs; j < sp_path.size(); ++j) {
        double n_spin = ((double) mults[j + pairs] - 1.0) / 2.0;
        double n_proj = (double) block_to_nzi(block, j + pairs) - n_spin;
        double delta_t_spin = (double) sp_path[j] - n_spin;
        c *= hashed_clebsh_gordan(t_spin, n_spin,t_spin + delta_t_spin,
                                  t_proj, n_proj);
        t_spin = t_spin + delta_t_spin;
        t_proj = t_proj + n_proj;

        if (c == 0.0) {
            return c;
        }
    }
    double coeff = sym2block[sym][0].coeff;
    return c / coeff;
}

unsigned int Indexes_P2::character_multiplication(unsigned int a, unsigned int b) {
    // A * A == A => (0 + 0) == 0 (mod 2)
    // A * B == B => (0 + 1) == 1 (mod 2)
    // B * B == A => (1 + 1) == 0 (mod 2)
    return (a + b) % 2;
}

void Indexes_P2::construct_branching_diagram() {
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
                unsigned int cum_symm = character_multiplication(symm_of_new_pair, bdi.representation);
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
        bd_from.push_back({1, sp2d.representation, v});
        tensor_product_spins(bd_from, bd_to, v.begin(), v.end());
        tensor_product_spins(bd_from, bd_to, mults.begin() + 2 * pairs, mults.end());
        bd_result.insert(bd_result.end(), bd_from.begin(), bd_from.end());
    }

    std::stable_sort(bd_result.begin(), bd_result.end());
    bds = std::move(bd_result);

    for (const Spin_P2_Diagramm & pmf : bds) {
        for (auto n : pmf.pairs_mults) {
            std::cout << n << " ";
        }
        std::cout << "| ";
        for (auto n : pmf.path) {
            std::cout << n << " ";
        }
        std::cout << "| representation : " << pmf.representation << " | back_mult: " << pmf.total_mult << std::endl;
    }
    std::cout << std::endl;

    construct_spin_boundaries();
}

int Indexes_P2::bds_total_mult(int i) const {
    return bds[i].total_mult;
}

unsigned int Indexes_P2::bds_total_repr(int i) const {
    return bds[i].representation;
}

unsigned long Indexes_P2::bds_size() const {
    return bds.size();
}

//        for (int mult : v) {
//            bd_to.clear();
//            for (const Spin_P2_Diagramm & bdi : bd_from) {
//                int t_mult = bdi.back_mult;
//                std::vector<unsigned int> t_full = bdi.path;
//                for (int k = std::abs(t_mult - mult) + 1; k < t_mult + mult; ++++k) {
//                    bd_to.push_back({k, bdi.representation, bdi.pairs_mults, t_full});
//                    bd_to.back().path.push_back((k - t_mult + mult) / 2);
//                }
//            }
//            std::swap(bd_from, bd_to);
//        }
//        for (unsigned int i = 2 * pairs; i < mults.size(); ++i) {
//            int mult = mults[i];
//            bd_to.clear();
//            for (const Spin_P2_Diagramm & bdi : bd_from) {
//                int t_mult = bdi.back_mult;
//                std::vector<unsigned int> t_full = bdi.path;
//                for (int k = std::abs(t_mult - mult) + 1; k < t_mult + mult; ++++k) {
//                    bd_to.push_back({k, bdi.representation, bdi.pairs_mults, t_full});
//                    bd_to.back().path.push_back((k - t_mult + mult) / 2);
//                }
//            }
//            std::swap(bd_from, bd_to);
//        }
