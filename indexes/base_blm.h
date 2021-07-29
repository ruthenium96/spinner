#pragma once
#include <vector>
#include <map>
#include <deque>
#include <algorithm>
#include <armadillo>
#include <boost/multi_array.hpp>

#include "quantum_numbers.h"
#include "../libs/wignerSymbols/include/wignerSymbols.h"


// TODO: для абелевых групп отлично сработает вместо path использовать
// TODO: условный вектор, размера N + (N - 1), где первые N -- это спины частиц
// TODO: остальные N-1 -- спины, получившиеся при N-1 сложении спинов
// TODO: соответственно, для каждой системы необходимо генерировать вектор инструкций,
// TODO: состоящий из пар (i, j), такой вектор на k-ом месте говорит,
// TODO: что в ячейку k+N вектора нужно записать результат сложения спинов, стоящих на местах i и j.

// TODO: продумать, как можно в схеме выше строить представления
// TODO: вектор, отмечающий блоки и их ориентацию (уместно для D2)?

class Sum_and_Lex {
public:
    unsigned long total_projection;
    unsigned long lex_index;

    bool operator<(const Sum_and_Lex &other) const {
        return this->total_projection < other.total_projection;
    }
};

// TODO: NaN coefficient for coeff = 0?
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
    // number of different representations in the group
    unsigned int num_of_repr;
    // size of group, i.e. number of different symmetry operations
    unsigned int group_size;

    int lex_to_nzi(unsigned long lex, int i) const;
    unsigned long nzi_to_lex(int nzi, int i) const;
    // it seems like there is no good way to build nzi => block, only nzs = > block
//    unsigned long nzi_to_block(int nzi, int i) const;
    unsigned long lex_to_block(unsigned long lex) const;
    unsigned long block_to_lex(unsigned long block) const;
    int block_to_nzi(unsigned long block, int i) const;
    std::vector<int> block_to_nzs(unsigned long block) const;
    unsigned long nzs_to_block(const std::vector<int>& nzs) const;
    unsigned long lex_ladder(unsigned long lex, int i, int ladder) const;
    unsigned long block_ladder(unsigned long block, int i, int ladder) const;

    arma::dmat construct_s2_transformation(int spin, unsigned int symmetry) const;

    const std::vector<Decomposition>& block_to_sym(unsigned long block) const;
    const std::vector<Decomposition>& sym_to_block(unsigned long sym) const;

    std::vector<std::vector<unsigned long>> sym_sum_boundaries;
    std::vector<std::vector<unsigned long>> sym_spin_boundaries;

    // symmetry-dependent function:
    virtual void construct_branching_diagram() = 0;

protected:
    mutable boost::multi_array<double, 5> CGs;
    double max_t_spin;
    double max_n_spin;

    virtual int bds_total_mult(int i) const = 0;
    virtual unsigned int bds_total_repr(int i) const = 0;
    // размер bds, то есть количество общих спинов, полученное при сложении всех спинов,
    // то есть количество неприводимых представлений, полученных при тензорном произведении начальных спинов из mults
    virtual unsigned long bds_size() const = 0;

    void resize_CGs();
    double hashed_clebsh_gordan(double l1, double l2, double l3, double m1, double m2) const;
    // symmetry-dependent function:
    virtual double total_coeff(unsigned long state, unsigned long sym) const = 0;

    std::vector<unsigned long> lex2block;
    std::vector<unsigned long> block2lex;
    std::vector<int> cumulative_product;

    void construct_cumulative_prod();

    std::vector<Sum_and_Lex> tensor_product();

    void construct_two_side_maps(const std::vector<Sum_and_Lex>& sum_and_lexes);

    //num_of_block_functions * num_of_repr (or less)
    std::vector<std::vector<Decomposition>> block2sym;
    //num_of_sym_functions * group_size (or any dividers of it)
    std::vector<std::vector<Decomposition>> sym2block;
    void construct_spin_boundaries();



    std::vector<std::array<int, 2>> spin_addition_scheme;
    std::deque<Quantum_Numbers> qns;
    void spin_addition();
    double total_coeff_new(unsigned long state, unsigned long sym) const;

    template <class SD, class IT>
    void tensor_product_spins(std::vector<SD> & bd_from, std::vector<SD> & bd_to,
                              const IT m_from, const IT m_to) {
        for (auto it_mult = m_from; it_mult != m_to; ++it_mult) {
            int mult = *it_mult;
            bd_to.clear();
            for (const SD & bdi : bd_from) {
                int t_mult = bdi.total_mult;
                for (int k = std::abs(t_mult - mult) + 1; k < t_mult + mult; ++++k) {
                    bd_to.push_back(bdi);
                    bd_to.back().total_mult = k;
                    bd_to.back().path.push_back((k - t_mult + mult) / 2);
                }
            }
            std::swap(bd_from, bd_to);
        }
    }
};
