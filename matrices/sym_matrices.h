#ifndef JULY_SYM_MATRICES_H
#define JULY_SYM_MATRICES_H

#include "matrices.h"

struct Symm_Matrices : Matrices {

    const Indexes_P2 & blm;

    void add_to_s_squared(unsigned int i, unsigned int j, double v) override {
        block_of_s_squared(i - current_start, j - current_start) += v;
    }

    Add_To_Matrix s_squared_ = static_cast<Add_To_Matrix>(&Symm_Matrices::add_to_s_squared);


    void add_to_hamiltonian(unsigned int i, unsigned int j, double v) override {
        block_of_hamiltonian(i - current_start, j - current_start) += v;
    }

    Add_To_Matrix hamiltonian_ = static_cast<Add_To_Matrix>(&Symm_Matrices::add_to_hamiltonian);


    Symm_Matrices(const Indexes_P2 &blm_,
                  const arma::dmat &js) : Matrices(blm_), js(js), blm(blm_) {
    }

    void construct_hamiltonian(unsigned int repr) {
        unsigned long current_length = current_end - current_start;
        block_of_hamiltonian.resize(current_length, current_length);
        block_of_hamiltonian.zeros(current_length, current_length);
        #pragma omp parallel for
        for (unsigned long sym = current_start; sym < current_end; ++sym) {
            for (int a = 0; a < blm.v_size; ++a) {
                for (int b = a + 1; b < blm.v_size; ++b) {
                    if (!std::isnan(js(a, b))) {
                        scalar_product_total(hamiltonian_, repr, -2 * js(a, b), sym, a, b);
                    }
                }
            }
        }
        std::cout << block_of_hamiltonian << std::endl;
    }


    void construct_s_squared(unsigned int repr) {
        unsigned long current_length = current_end - current_start;
        block_of_s_squared.resize(current_length, current_length);
        block_of_s_squared.zeros(current_length, current_length);
        double diagonal_sum = 0;
        for (double s : blm.spins) {
            diagonal_sum += s * (s + 1);
        }

#pragma omp parallel for
        for (unsigned long sym = current_start; sym < current_end; ++sym) {
            add_to_s_squared(sym, sym, diagonal_sum);
            for (int a = 0; a < blm.v_size; ++a) {
                for (int b = a + 1; b < blm.v_size; ++b) {
                    if (a != b) {
                        scalar_product_total(s_squared_, repr, 2, sym, a, b);
                    }
                }
            }
        }
        std::cout << block_of_s_squared << std::endl;
    }


    void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                            arma::vec & degeneracy) override {
        eigval.reset();
        s_squared_new_basis_vector.reset();
        degeneracy.reset();
        for (unsigned int r = 0; r < blm.num_of_repr; ++r) {
            current_start = blm.sym_sum_boundaries[r].front();
            current_end = blm.sym_sum_boundaries[r].back();

            double a = 1/sqrt(2);
            double b = -a;

            construct_hamiltonian(r);

//            arma::dmat rotation_one;
//            rotation_one.resize(16, 16);
//            rotation_one(0, 0) = 1;
//            rotation_one(1, 1) = 1;
//            rotation_one(2, 2) = 1;
//            rotation_one(3, 3) = a; rotation_one(3, 12) = a;
//            rotation_one(4, 3) = b; rotation_one(4, 12) = a;
//            rotation_one(5, 4) = 1;
//            rotation_one(6, 5) = a; rotation_one(6, 13) = a;
//            rotation_one(8, 5) = b; rotation_one(8, 13) = a;
//            rotation_one(7, 6) = a; rotation_one(7, 14) = a;
//            rotation_one(9, 6) = b; rotation_one(9, 14) = a;
//            rotation_one(10, 7) = 1;
//            rotation_one(11, 8) = a; rotation_one(11, 15) = a;
//            rotation_one(12, 8) = b; rotation_one(12, 15) = a;
//            rotation_one(13, 9) = 1;
//            rotation_one(14, 10) = 1;
//            rotation_one(15, 11) = 1;

//            std::cout << rotation_one << std::endl;
//
//            std::cout << rotation_one.t() * block_of_hamiltonian * rotation_one << std::endl;





            arma::mat eigvec;
            arma::vec new_eigval;
            arma::eig_sym(new_eigval, eigvec, block_of_hamiltonian);
            eigval = arma::join_cols(eigval, new_eigval);
            block_of_hamiltonian.reset();
//
//            construct_s_squared(r);
//            arma::mat s_squared_new_basis = eigvec.t() * block_of_s_squared * eigvec;
//            arma::vec new_s_squared_new_basis_vector = s_squared_new_basis.diag();
//            s_squared_new_basis_vector = arma::join_cols(s_squared_new_basis_vector, new_s_squared_new_basis_vector);
//            degeneracy = arma::ones(size(eigval));
        }
        eigval -= eigval.min();
        std::cout << eigval << std::endl;
    }

    // "factor" is the factor before \sum_{i<j}(S_i; S_j)
    // for HDvV factor should be -2J
    // for S^2 factor should be 2:
    // S_{T}^2 = (\sum_{i}S_i)^2 = \sum_{i}S_i^2 + 2 * \sum_{i<j}(S_i; S_j)
//    void scalar_product_total_sym(Add_To_Matrix add_to_matrix, unsigned int r, const double factor,
//                                  const unsigned long sym, const int a, const int b) {
//        for (const Decomposition & d : blm.sym2block[sym]) {
//            int na = blm.block_to_nzi(d.index, a);
//            int nb = blm.block_to_nzi(d.index, b);
//
//            double decomposition_factor = d.coeff * blm.block2sym[d.index][r].coeff;
//
//            // Saz Sbz
//            (this->*add_to_matrix)(sym, sym, (na - blm.mults[a]) * (nb - blm.mults[b]) * factor * decomposition_factor);
//            // Sa+ Sb-
//            scalar_product_plus_minus_sym(add_to_matrix, r, factor, sym, d, a, b, na, nb);
//            // Sa- Sb+
//            scalar_product_plus_minus_sym(add_to_matrix, r, factor, sym, d, b, a, nb, na);
//        }
//    }

//    void scalar_product_plus_minus_sym(Add_To_Matrix add_to_matrix, unsigned int r, const double factor,
//                                       const unsigned long sym, const Decomposition & d,
//                                       const int a, const int b,
//                                       const int na, const int nb) {
//        if (na == blm.mults[a] - 1 || nb == 0) {
//            return;
//        }
//        unsigned long new_block = blm.block_ladder(blm.block_ladder(d.index, a, +1), b, -1);
//
//        if (r < blm.block2sym[new_block].size()) {
//            unsigned long new_sym = blm.block2sym[new_block][r].index;
//            double decomposition_factor = d.coeff * blm.block2sym[new_block][r].coeff;
//
//            // projection m = number n - spin S
//            // so S(S+1)-m(m+1) = (2S-n)(n+1)
//            // so S(S+1)-m(m-1) = n(2S+1-n)
//            double factor_a = (2 * blm.mults[a] - na) * (na + 1);
//            double factor_b = nb * (2 * blm.mults[b] + 1 - nb);
//
//            (this->*add_to_matrix)(sym, new_sym, 0.5 * sqrt(factor_a * factor_b) * factor * decomposition_factor);
//        }
//    }


    const arma::dmat &js;
    arma::dmat block_of_hamiltonian;
    arma::dmat block_of_s_squared;
    unsigned int current_start;
    unsigned int current_end;
};

#endif
