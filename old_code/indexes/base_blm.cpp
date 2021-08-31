#include "base_blm.h"

#include <iostream>
#include <utility>

Indexes::Indexes(std::vector<int> mults_): mults(std::move(mults_)) {

    spins.resize(mults.size());
    for (int i = 0; i < mults.size(); ++i) {
        spins[i] = (mults[i] - 1) / 2.0;
    }

    // Constructing of lex <-> block transformation
    std::vector<Sum_and_Lex> sum_and_lexes = tensor_product();
    std::stable_sort(sum_and_lexes.begin(), sum_and_lexes.end());
    construct_two_side_maps(sum_and_lexes);
}

int Indexes::lex_to_nzi(unsigned long lex, int i) const {
    return (lex % cumulative_product[i]) / cumulative_product[i + 1];
}

unsigned long Indexes::nzi_to_lex(int nzi, int i) const{
    return nzi * cumulative_product[i + 1];
}

unsigned long Indexes::lex_to_block(unsigned long lex) const{
    return lex2block[lex];
}

unsigned long Indexes::block_to_lex(unsigned long block) const{
    return block2lex[block];
}

int Indexes::block_to_nzi(unsigned long block, int i) const{
    return lex_to_nzi(block_to_lex(block), i);
}

std::vector<int> Indexes::block_to_nzs(unsigned long block) const{
    std::vector<int> nzs(v_size);
    for (int i = 0; i < v_size; ++i) {
        nzs[i] = lex_to_nzi(block_to_lex(block), i);
    }
    return std::move(nzs);
}

unsigned long Indexes::nzs_to_block(const std::vector<int>& nzs) const{
    unsigned long lex = 0;
    for (int i = 0; i < v_size; ++i) {
        lex += nzi_to_lex(nzs[i], i);
    }
    return lex2block[lex];
}

void Indexes::construct_cumulative_prod() {
    cumulative_product.resize(v_size + 1, 0);
    cumulative_product[0] = tensor_size;
    for (int k = 1; k < v_size + 1; ++k) {
        cumulative_product[k] = cumulative_product[k - 1] / mults[k - 1];
    }
}

// NB: one should check projection value before ladding
unsigned long Indexes::lex_ladder(unsigned long lex, int i, int ladder) const{
    return lex + ladder * cumulative_product[i + 1];
}

// NB: one should check projection value before ladding
unsigned long Indexes::block_ladder(unsigned long block, int i, int ladder) const{
    return lex_to_block(lex_ladder(block_to_lex(block), i, ladder));
}

std::vector<Sum_and_Lex> Indexes::tensor_product() {
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

void Indexes::construct_two_side_maps(const std::vector<Sum_and_Lex>& sum_and_lexes) {
    lex2block.resize(tensor_size, 0);
    block2lex.resize(tensor_size, 0);

    unsigned long current_sum = 0;

    for (unsigned long block = 0; block < tensor_size; ++block) {
        unsigned long lex = sum_and_lexes[block].lex_index;
        block2lex[block] = lex;
        lex2block[lex] = block;
        current_sum = std::max(current_sum, sum_and_lexes[block].total_projection);
    }
}

arma::dmat Indexes::construct_s2_transformation(int spin, unsigned int repr) const {
    // boundaries in state-basis
    unsigned long br_start = sym_spin_boundaries[repr][spin];
    unsigned long br_end = sym_spin_boundaries[repr][spin + 1];
    unsigned long br_len = br_end - br_start;

    // boundaries in projection-basis
    unsigned long current_start = sym_sum_boundaries[repr][spin];
    unsigned long current_end = sym_sum_boundaries[repr][spin + 1];
    unsigned long bl_len = current_end - current_start;

    arma::dmat trans(br_len, bl_len);
    #pragma omp parallel for
    for (unsigned long state = br_start; state < br_end; ++state) {
        for (unsigned long sym = current_start; sym < current_end; ++sym) {
            trans(state - br_start, sym - current_start) = total_coeff_new(state, sym);
        }
    }

//    std::cout << trans << std::endl;

    return std::move(trans);
}

// Clebsch-Gordan function has six parameters. It's hard to hash it.
// But the sixth parameter (projection of total spin) should be equal the sum of other
// projections, otherwise coefficient is zero. So we can hash only five parameters.
// Total spin, by a triangle law, possess values from |S1 - S2| to (S1 + S2),
// this fact also can be used for optimization.
void Indexes::resize_CGs() {
    boost::multi_array<double, 5>::extent_gen extents;
    int max_t_spin_number = (int) (2 * max_t_spin + 1);
    int max_n_spin_number = (int) (2 * max_n_spin + 1);
    // a triangle law:
    // TODO: test it and maybe refactor!
    int max_s_spin_number = (int) (4 * max_n_spin + 1);
//    CGs.resize(extents[max_t_spin_number][max_n_spin_number][max_t_spin_number][2 * max_t_spin_number + 1][2 * max_n_spin_number + 1]);
    CGs.resize(extents[max_t_spin_number][max_n_spin_number][max_s_spin_number][2 * max_t_spin_number + 1][2 * max_n_spin_number + 1]);
    std::fill_n(CGs.data(), CGs.num_elements(), NAN);
}

double Indexes::hashed_clebsh_gordan(double l1, double l2, double l3,
                                     double m1, double m2) const {
    // a triangle law:
    if (l3 > l1 + l2 || l3 < std::abs(l1 - l2)) {
        return 0;
    }
    int il1 = (int) (2 * l1);
    int il2 = (int) (2 * l2);

    // a triangle law:
    // TODO: test it and maybe refactor!
    int min_l = std::abs(il1 - il2);
    int il3 = (int) (2 * l3) - min_l;

    int im1 = (int) (2 * (m1 + max_t_spin));
    int im2 = (int) (2 * (m2 + max_n_spin));
    double& curr_value = CGs[il1][il2][il3][im1][im2];
    if (std::isnan(curr_value)) {
        curr_value = WignerSymbols::clebschGordan(l1, l2, l3,
                                                  m1, m2,m1 + m2);
    }
    return curr_value;

}

// Construct sym_spin_boundaries vector.
void Indexes::construct_spin_boundaries() {

    // TODO: аналогичная переменная есть для CGs.
    int max_mult = (2.0 * std::accumulate(spins.begin(), spins.end(), 0.0)) + 1;

    int current_mult = max_mult;
    unsigned int current_repr = 0;

    // number of different values of possible spin:
    int diff_values_of_spin = ceil((double) max_mult / 2.0);

//    sym_spin_boundaries.resize(num_of_repr, std::vector<unsigned long>(diff_values_of_spin + 1, bds_size()));
    sym_spin_boundaries.resize(num_of_repr, std::vector<unsigned long>(diff_values_of_spin + 1, qns.size()));

    for (unsigned long i = 0; i < qns.size(); ++i) {
//        if (bds_total_repr(i) != current_repr) {
        if (qns[i].representation != current_repr) {

            for (; current_mult > -2; ----current_mult) {
                int index = (max_mult - current_mult) / 2;
                sym_spin_boundaries[current_repr][index] = i;
            }
            current_mult = max_mult;
            ++current_repr;
        }
//        for (; (bds_total_mult(i) <= current_mult) && (current_mult > 0); ----current_mult) {
        for (; (qns[i].back_mult() <= current_mult) && (current_mult > 0); ----current_mult) {

            int index = (max_mult - current_mult) / 2;
            sym_spin_boundaries[current_repr][index] = i;
        }
    }
        std::cout << "sym_sum_boundaries" << std::endl;

        for (const auto& v : sym_sum_boundaries) {
            for (auto n : v) {
                std::cout << n << " ";
            }
            std::cout << std::endl;
        }

    std::cout << "sym_spin_boundaries" << std::endl;

    for (const auto& v : sym_spin_boundaries) {
        for (auto n : v) {
            std::cout << n << " ";
        }
        std::cout << std::endl;
    }
}

const std::vector<Decomposition> &Indexes::block_to_sym(unsigned long block) const {
    return block2sym[block];
}

const std::vector<Decomposition> &Indexes::sym_to_block(unsigned long sym) const {
    return sym2block[sym];
}

void Indexes::spin_addition() {
    qns.emplace_back(mults, 0);

    int step = 0;
    std::array<int, 2> pair = spin_addition_scheme[step];
    while (step < v_size - 1) {
        int mult_one = qns.front()(pair[0]);
        int mult_two = qns.front()(pair[1]);
        for (int mult_new = std::abs(mult_one - mult_two) + 1; mult_new < mult_one + mult_two; ++++mult_new) {
            qns.push_back(qns.front());
            qns.back().push_back(mult_new);
        }

        qns.pop_front();
        if (qns.front().size() > step) {
            ++step;
            pair = spin_addition_scheme[step];
        }
    }
    std::stable_sort(qns.begin(), qns.end());
//    for (auto v : qns) {
//        v.print();
//    }
}

double Indexes::total_coeff_new(unsigned long state, unsigned long sym) const {
    // TODO: It works only for P1 symmetry now.
    unsigned long block = sym;
    const Quantum_Numbers & qni = qns[state];
    double c = 1;
    // TODO: temporary (and mutable) vector<vector> or arma::dmat for cum_proj
    std::vector<double> cum_proj(qni.size(), NAN);
    for (int i = 0; i < qni.size(); ++i) {
        int n_of_spin_one = spin_addition_scheme[i][0];
        int n_of_spin_two = spin_addition_scheme[i][1];
        // TODO: implement spin_to_mult & mult_to_spin?
        double spin_one = ((double) qni(n_of_spin_one) - 1.0) / 2.0;
        double spin_two = ((double) qni(n_of_spin_two) - 1.0) / 2.0;
        double spin_sum = ((double) qni(v_size + i) - 1.0) / 2.0;
        double proj_one, proj_two;
        // TODO: абстрагировать одинаковый код
        if (n_of_spin_one < v_size) {
            proj_one = (double) block_to_nzi(block, n_of_spin_one) - spin_one;
        } else {
            proj_one = cum_proj[n_of_spin_one - v_size];
        }
        if (n_of_spin_two < v_size) {
            proj_two = (double) block_to_nzi(block, n_of_spin_two) - spin_two;
        } else {
            proj_two = cum_proj[n_of_spin_two - v_size];
        }
        cum_proj[i] = proj_one + proj_two;
        c *= hashed_clebsh_gordan(spin_one, spin_two, spin_sum,
                                  proj_one, proj_two);
        if (c == 0.0) {
            return c;
        }
    }
    return c;
}

