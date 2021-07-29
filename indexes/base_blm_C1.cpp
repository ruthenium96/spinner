#include "base_blm_C1.h"

Indexes_P1::Indexes_P1(std::vector<int> mults_) : Indexes(mults_) {
    num_of_repr = 1;
    group_size = 1;

    dim_of_repr.resize(num_of_repr);
    dim_of_repr[0] = tensor_size;
    sym_sum_boundaries.resize(num_of_repr);

    block2sym.resize(tensor_size);
    sym2block.resize(tensor_size);
    construct_symm_decomposition();

    max_t_spin = std::accumulate(spins.begin(), spins.end(), 0.0);
    max_n_spin = (*std::max_element(spins.begin(), spins.end()));
    resize_CGs();
}

int Indexes_P1::bds_total_mult(int i) const {
    return bds[i].total_mult;
}

unsigned int Indexes_P1::bds_total_repr(int i) const {
    return 0;
}

unsigned long Indexes_P1::bds_size() const {
    return bds.size();
}

double Indexes_P1::total_coeff(unsigned long state, unsigned long sym) const {
    unsigned long block = sym;
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
        t_spin = t_spin + delta_t_spin;
        t_proj = t_proj + n_proj;

        if (c == 0.0) {
            return c;
        }
    }
    return c;
}

void Indexes_P1::construct_branching_diagram() {

    spin_addition_scheme.resize(v_size - 1);

    spin_addition_scheme[0][0] = 0;
    spin_addition_scheme[0][1] = 1;
    for (int i = 2; i < v_size; ++i) {
        spin_addition_scheme[i - 1][0] = (v_size - 2) + i;
        spin_addition_scheme[i - 1][1] = i;
    }

//    for (auto v : spin_addition_scheme) {
//        std::cout << v[0] << " " << v[1] << std::endl;
//    }

    spin_addition();

    // TODO: rewrite without pushbacks!
    std::vector<Spin_P1_Diagramm> bd_from;
    std::vector<Spin_P1_Diagramm> bd_to;
    bd_from.push_back({1});

    tensor_product_spins(bd_from, bd_to, mults.begin(), mults.end());

    std::stable_sort(bd_from.begin(), bd_from.end());
    bds = std::move(bd_from);

    construct_spin_boundaries();
}

void Indexes_P1::construct_symm_decomposition() {
    std::vector<unsigned int> counters;
    counters.resize(num_of_repr);
    counters[0] = 0;

    unsigned long current_sum = UINT64_MAX;
    for (unsigned long block = 0; block < tensor_size; ++block) {
        std::vector<int> nzs = block_to_nzs(block);
        unsigned long sum = std::accumulate(nzs.begin(), nzs.end(), 0);
        if (sum != current_sum) {
            sym_sum_boundaries[0].push_back(counters[0]);
            current_sum = sum;
        }
        if (block2sym[block].empty()) {
            block2sym[block] = std::vector<Decomposition>(1);
            block2sym[block][0] = {counters[0], 1.0};

            sym2block[counters[0]] = std::vector<Decomposition>(1);
            sym2block[counters[0]][0] = {block, 1.0};

            ++counters[0];
        }
    }
    sym_sum_boundaries[0].push_back(tensor_size);
}

//    for (int mult : mults) {
//        bd_to.clear();
//        for (const Spin_P1_Diagramm & bdi : bd_from) {
//            int t_mult = bdi.back_mult;
//            std::vector<unsigned int> t_full = bdi.path;
//            for (int k = std::abs(t_mult - mult) + 1; k < t_mult + mult; ++++k) {
//                bd_to.push_back({k, t_full});
//                bd_to.back().path.push_back((k - t_mult + mult) / 2);
//            }
//        }
//        std::swap(bd_from, bd_to);
//    }
