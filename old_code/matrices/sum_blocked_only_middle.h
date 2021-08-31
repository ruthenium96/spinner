#ifndef JULY_SUM_BLOCKED_ONLY_MIDDLE_H
#define JULY_SUM_BLOCKED_ONLY_MIDDLE_H

#include "base_matrices.h"

// build and diagonalize central block of Hamiltonian.
// then multiply formulas to multiplicity
struct Proj_Middle_Blocked_Matrices : Matrices {
    const Indexes_P2 & blm;
    void add_to_s_squared(unsigned int i, unsigned int j, double v) override {
        block_of_s_squared(i - current_start, j - current_start) += v;
    }

    Add_To_Matrix s_squared_ = static_cast<Add_To_Matrix>(&Proj_Middle_Blocked_Matrices::add_to_s_squared);


    void add_to_hamiltonian(unsigned int i, unsigned int j, double v) override {
        block_of_hamiltonian(i - current_start, j - current_start) += v;
    }

    Add_To_Matrix hamiltonian_ = static_cast<Add_To_Matrix>(&Proj_Middle_Blocked_Matrices::add_to_hamiltonian);

    Proj_Middle_Blocked_Matrices(const Indexes_P2 &blm_,
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
//        std::cout << block_of_hamiltonian << std::endl;
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
//        std::cout << block_of_s_squared << std::endl;
    }

    void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                            arma::vec & degeneracy) override {
        eigval.reset();
        s_squared_new_basis_vector.reset();
        degeneracy.reset();
        for (unsigned int r = 0; r < blm.num_of_repr; ++r) {
            int middle = blm.sym_sum_boundaries[r].size() / 2 - 1;
            current_start = blm.sym_sum_boundaries[r][middle];
            current_end = blm.sym_sum_boundaries[r][middle + 1];

            construct_hamiltonian(r);
            arma::mat eigvec;
            arma::vec new_eigval;
            arma::eig_sym(new_eigval, eigvec, block_of_hamiltonian);
            eigval = arma::join_cols(eigval, new_eigval);
            block_of_hamiltonian.reset();

            construct_s_squared(r);
            arma::mat s_squared_new_basis = eigvec.t() * block_of_s_squared * eigvec;
            arma::vec new_s_squared_new_basis_vector = s_squared_new_basis.diag();

            arma::vec new_degeneracy = sqrt(4 * new_s_squared_new_basis_vector + arma::ones(size(new_s_squared_new_basis_vector)));

            s_squared_new_basis_vector = arma::join_cols(s_squared_new_basis_vector, new_s_squared_new_basis_vector);
            degeneracy = arma::join_cols(degeneracy, new_degeneracy);

        }
        eigval -= eigval.min();
//        std::cout << eigval << std::endl << std::endl << s_squared_new_basis_vector << std::endl << std::endl << degeneracy << std::endl;
    }

    const arma::dmat &js;
    arma::dmat block_of_hamiltonian;
    arma::dmat block_of_s_squared;
    unsigned int current_start;
    unsigned int current_end;
};

#endif //JULY_SUM_BLOCKED_ONLY_MIDDLE_H