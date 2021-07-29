#ifndef JULY_SUM_S2_BLOCKED_H
#define JULY_SUM_S2_BLOCKED_H

#include <algorithm>
#include "base_matrices.h"

struct S2_Proj_Blocked_Matrices : Matrices {
    const Indexes & blm;
    void add_to_s_squared(unsigned int i, unsigned int j, double v) override {
    }
    Add_To_Matrix s_squared_ = static_cast<Add_To_Matrix>(&S2_Proj_Blocked_Matrices::add_to_s_squared);

    void add_to_hamiltonian(unsigned int i, unsigned int j, double v) override {
        block_of_hamiltonian(i - current_start, j - current_start) += v;
    }
    Add_To_Matrix hamiltonian_ = static_cast<Add_To_Matrix>(&S2_Proj_Blocked_Matrices::add_to_hamiltonian);

    S2_Proj_Blocked_Matrices(const Indexes &blm_,
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

    void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                            arma::vec & degeneracy) override {
        eigval.reset();
        s_squared_new_basis_vector.reset();
        degeneracy.reset();
        for (unsigned int r = 0; r < blm.num_of_repr; ++r) {
            for (int s = 0; s < blm.sym_spin_boundaries[r].size() - 1; ++s) {
                // TOOD: implement way
                double multiplicity = (2.0 * (std::accumulate(blm.spins.begin(), blm.spins.end(), 0.0) - s)) + 1;
                double s_squared = (multiplicity * multiplicity - 1.0) / 4.0;

                current_start = blm.sym_sum_boundaries[r][s];
                current_end = blm.sym_sum_boundaries[r][s + 1];

                construct_hamiltonian(r);

                arma::dmat trans = blm.construct_s2_transformation(s, r);

//                std::cout << "s = " << s << " r = " << r << std::endl << trans << std::endl;

                arma::dmat new_hamiltonian = trans * block_of_hamiltonian * trans.t();

                arma::vec new_eigval;
                arma::eig_sym(new_eigval, new_hamiltonian);
                eigval = arma::join_cols(eigval, new_eigval);

//                std::cout << block_of_hamiltonian << std::endl;
//                std::cout << new_eigval << std::endl;

                block_of_hamiltonian.reset();

                degeneracy = arma::join_cols(degeneracy, multiplicity * arma::ones(new_eigval.n_rows, 1));
                s_squared_new_basis_vector = arma::join_cols(s_squared_new_basis_vector, s_squared * arma::ones(new_eigval.n_rows, 1));
            }
        }
        eigval -= eigval.min();
//        std::cout << eigval << std::endl << std::endl << s_squared_new_basis_vector << std::endl << std::endl << degeneracy << std::endl;
    }

    const arma::dmat &js;
    arma::dmat block_of_hamiltonian;
    unsigned int current_start;
    unsigned int current_end;
};


#endif //JULY_SUM_S2_BLOCKED_H