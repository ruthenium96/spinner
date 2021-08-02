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


// It is a "branching diagram", i.e. the diagram produced by step-by-step addition of mults.
// 1) S_1
// 2) (S_1) + S_2
// 3) (S_1 + S_2) + S_3
// i) (S_1 + S_2 + ... + S_{i-1}) + S_i
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