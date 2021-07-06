#include "base_blm.h"

#include <iostream>
#include <utility>

Block_Lex_Map::Block_Lex_Map(std::vector<int> mults_): mults(std::move(mults_)) {

    spins.resize(mults.size());
    for (int i = 0; i < mults.size(); ++i) {
        spins[i] = (mults[i] - 1) / 2.0;
    }

    std::vector<Sum_and_Lex> sum_and_lexes = tensor_product();
    std::stable_sort(sum_and_lexes.begin(), sum_and_lexes.end());
    construct_two_side_maps(sum_and_lexes);
}

int Block_Lex_Map::lex_to_nzi(unsigned long lex, int i) const {
    return (lex % cumulative_product[i]) / cumulative_product[i + 1];
}

unsigned long Block_Lex_Map::nzi_to_lex(int nzi, int i) const{
    return nzi * cumulative_product[i + 1];
}

unsigned long Block_Lex_Map::nzi_to_block(int nzi, int i) const{
    return lex2block[nzi * cumulative_product[i + 1]];
}

unsigned long Block_Lex_Map::lex_to_block(unsigned long lex) const{
    return lex2block[lex];
}

unsigned long Block_Lex_Map::block_to_lex(unsigned long block) const{
    return block2lex[block];
}

int Block_Lex_Map::block_to_nzi(unsigned long block, int i) const{
    return lex_to_nzi(block_to_lex(block), i);
}

std::vector<int> Block_Lex_Map::block_to_nzs(unsigned long block) const{
    std::vector<int> nzs(v_size);
    for (int i = 0; i < v_size; ++i) {
        nzs[i] = lex_to_nzi(block_to_lex(block), i);
    }
    return std::move(nzs);
}

unsigned long Block_Lex_Map::nzs_to_block(const std::vector<int>& nzs) const{
    unsigned long lex = 0;
    for (int i = 0; i < v_size; ++i) {
        lex += nzi_to_lex(nzs[i], i);
    }
    return lex2block[lex];
}

void Block_Lex_Map::construct_cumulative_prod() {
    cumulative_product.resize(v_size + 1, 0);
    cumulative_product[0] = tensor_size;
    for (int k = 1; k < v_size + 1; ++k) {
        cumulative_product[k] = cumulative_product[k - 1] / mults[k - 1];
    }
}

// one should check projection value before ladding
unsigned long Block_Lex_Map::lex_ladder(unsigned long lex, int i, int ladder) const{
    return lex + ladder * cumulative_product[i + 1];
}

unsigned long Block_Lex_Map::block_ladder(unsigned long block, int i, int ladder) const{
    return lex_to_block(lex_ladder(block_to_lex(block), i, ladder));
}

std::vector<Sum_and_Lex> Block_Lex_Map::tensor_product() {
    tensor_size = 1;
    for (int mult : mults) {
        tensor_size *= mult;
    }
    v_size = mults.size();

    construct_cumulative_prod();
    std::vector<Sum_and_Lex> sum_and_lexes(tensor_size);

    #pragma omp parallel for
    for (unsigned long j = 0; j < tensor_size; ++j) {
        unsigned long local_sum = 0;
        for (int i = 0; i < v_size; ++i) {
            local_sum += lex_to_nzi(j, i);
        }
        sum_and_lexes[j] = {local_sum, j};
    }
    return std::move(sum_and_lexes);
}

void Block_Lex_Map::construct_two_side_maps(const std::vector<Sum_and_Lex>& sum_and_lexes) {
    lex2block.resize(tensor_size, 0);
    block2lex.resize(tensor_size, 0);
    // TODO: rewrite block_boundaries part without push_back's
    // TODO: отделить построение block_boundaries от остального, это чуть-чуть всё замедлит, но зато удачно разделит
    // TODO: заодно добавить бы проверку того, что block_boundaries построен

    unsigned long current_sum = 0;
    block_boundaries.push_back(current_sum);

    for (unsigned long block = 0; block < tensor_size; ++block) {
        unsigned long lex = sum_and_lexes[block].lex_index;
        block2lex[block] = lex;
        // memory random access:
        lex2block[lex] = block;
        if (sum_and_lexes[block].sum > current_sum) {
            current_sum = sum_and_lexes[block].sum;
            block_boundaries.push_back(block);
        }
    }
    block_boundaries.push_back(tensor_size);
}


