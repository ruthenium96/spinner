#include "base_blm_C2.h"

Indexes_P2::Indexes_P2(std::vector<int> mults_,
                       unsigned int pairs_) : Indexes(std::move(mults_)), pairs(pairs_) {
    dim_of_repr.resize(num_of_repr);
    counters.resize(num_of_repr);
    sym_sum_boundaries.resize(num_of_repr);
    block2sym.resize(tensor_size);
    sym2block.resize(tensor_size);
    calculate_dim_of_repr();
    counters[0] = 0;
    counters[1] = dim_of_repr[0];
    contruct_symm_decomposition();

    for (const auto& vd : block2sym) {
        for (auto d : vd) {
            std::cout << "(" << d.coeff << " * |" << d.index << ">) + ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for (const auto& vd : sym2block) {
        for (auto d : vd) {
            std::cout << "(" << d.coeff << " * |" << d.index << ">) + ";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;

    for (const auto& repr : sym_sum_boundaries) {
        for (auto n : repr) {
            std::cout << n << " ";
        }
        std::cout << std::endl;
    }
}

void Indexes_P2::contruct_symm_decomposition() {
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
//    for (auto nz : nzs) {
//        std::cout << nz << " ";
//    }
//    std::cout << std::endl;
    for (int k = 0; k < pairs; ++k) {
        std::swap(nzs[2 * k], nzs[2 * k + 1]);
    }
//    for (auto nz : nzs) {
//        std::cout << nz << " ";
//    }
//    std::cout << std::endl << std::endl;
    return nzs_to_block(nzs);
}
