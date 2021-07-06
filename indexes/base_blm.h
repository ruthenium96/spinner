#pragma once
#include <vector>
#include <map>
#include <algorithm>
#include "armadillo"

class Sum_and_Lex {
public:
    unsigned long sum;
    unsigned long lex_index;

    bool operator<(const Sum_and_Lex &other) const {
        return this->sum < other.sum;
    }
};

struct Block_Lex_Map {
public:
    explicit Block_Lex_Map(std::vector<int> mults_);
    Block_Lex_Map(const Block_Lex_Map&) = delete;
    Block_Lex_Map& operator= (const Block_Lex_Map&) = delete;
    Block_Lex_Map(Block_Lex_Map&&) noexcept = delete;
    Block_Lex_Map& operator= (Block_Lex_Map&&) = delete;
    ~Block_Lex_Map() noexcept = default;

    int tensor_size;
    unsigned long v_size;
    const std::vector<int> mults;
    std::vector<double> spins;

    // TODO: функция-обёртка для только константного доступа?
    std::vector<unsigned long> block_boundaries;

    int lex_to_nzi(unsigned long lex, int i) const;
    unsigned long nzi_to_lex(int nzi, int i) const;
    unsigned long nzi_to_block(int nzi, int i) const;
    unsigned long lex_to_block(unsigned long lex) const;
    unsigned long block_to_lex(unsigned long block) const;
    int block_to_nzi(unsigned long block, int i) const;
    std::vector<int> block_to_nzs(unsigned long block) const;
    unsigned long nzs_to_block(const std::vector<int>& nzs) const;
    unsigned long lex_ladder(unsigned long lex, int i, int ladder) const;
    unsigned long block_ladder(unsigned long block, int i, int ladder) const;

private:
    std::vector<unsigned long> lex2block;
    std::vector<unsigned long> block2lex;
    std::vector<int> cumulative_product;
    std::vector<std::array<int, 2>> vec_of_arr;

    void construct_cumulative_prod();

    std::vector<Sum_and_Lex> tensor_product();

    void construct_two_side_maps(const std::vector<Sum_and_Lex>& sum_and_lexes);
};
