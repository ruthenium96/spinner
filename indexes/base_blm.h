#pragma once
#include <vector>
#include <map>
#include <algorithm>
#include "armadillo"

class Sum_and_Lex {
public:
    unsigned long total_projection;
    unsigned long lex_index;

    bool operator<(const Sum_and_Lex &other) const {
        return this->total_projection < other.total_projection;
    }
};

struct Decomposition {
    unsigned long index;
    double coeff;
};


struct Indexes {
public:

    explicit Indexes(std::vector<int> mults_);
    Indexes(const Indexes&) = delete;
    Indexes& operator= (const Indexes&) = delete;
    Indexes(Indexes&&) noexcept = delete;
    Indexes& operator= (Indexes&&) = delete;
    ~Indexes() noexcept = default;

    int tensor_size;
    unsigned long v_size;
    const std::vector<int> mults;
    std::vector<double> spins;

    // TODO: функция-обёртка для только константного доступа?
//    std::vector<unsigned long> block_boundaries;

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

    std::vector<std::vector<Decomposition>> block2sym;
    //num_of_block_functions * num_of_repr (or less)

    std::vector<std::vector<Decomposition>> sym2block;
    //num_of_sym_functions * group_size (or any dividers of it)


private:
    std::vector<unsigned long> lex2block;
    std::vector<unsigned long> block2lex;
    std::vector<int> cumulative_product;

    void construct_cumulative_prod();

    std::vector<Sum_and_Lex> tensor_product();

    void construct_two_side_maps(const std::vector<Sum_and_Lex>& sum_and_lexes);
};
