#ifndef JULY_SUM_S_SQUARED_BLOCKED_H
#define JULY_SUM_S_SQUARED_BLOCKED_H

#include <algorithm>
#include "base_matrices.h"

struct Only_J_S_Squared : Matrices {
public:
    arma::dmat block_of_hamiltonian;

    void add_to_s_squared(unsigned int i, unsigned int j, double v) override {

    }

    void add_to_hamiltonian(unsigned int i, unsigned int j, double v) override {
        block_of_hamiltonian(i - current_start, j - current_start) += v;
    }

    Add_To_Matrix hamiltonian_ = static_cast<Add_To_Matrix>(&Only_J_S_Squared::add_to_hamiltonian);

    Only_J_S_Squared(const Indexes &blm,
                     const arma::dmat &js) : Matrices(blm), js(js) {}

    void construct_hamiltonian() {
        unsigned long current_length = current_end - current_start;
        block_of_hamiltonian.resize(current_length, current_length);
        block_of_hamiltonian.zeros(current_length, current_length);
        #pragma omp parallel for
        for (unsigned long block = current_start; block < current_end; ++block) {
            for (int a = 0; a < blm.v_size; ++a) {
                for (int b = a + 1; b < blm.v_size; ++b) {
                    if (!std::isnan(js(a, b))) {
                        scalar_product_total(hamiltonian_, -2 * js(a, b), block, a, b);
                    }
                }
            }
        }
    }

    void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                            arma::vec & degeneracy) override{
        eigval.reset();
        s_squared_new_basis_vector.reset();
        degeneracy.reset();
        for (int i = 0; i < blm.branch_boundaries.size() - 1; ++i) {
            double multiplicity = blm.bds[blm.branch_boundaries[i]].total_mult;
            double s_squared = (multiplicity * multiplicity - 1.0) / 4.0;
            current_start = blm.block_boundaries[i];
            current_end = blm.block_boundaries[i + 1];
            arma::dmat trans = blm.construct_transformation(i);
            construct_hamiltonian();
            arma::dmat blocked_hamiltonian = trans * block_of_hamiltonian * trans.t();

            arma::vec new_eigval;
            arma::eig_sym(new_eigval, blocked_hamiltonian);
            eigval = arma::join_cols(eigval, new_eigval);
            degeneracy = arma::join_cols(degeneracy, multiplicity * arma::ones(new_eigval.n_rows, 1));
            s_squared_new_basis_vector = arma::join_cols(s_squared_new_basis_vector, s_squared * arma::ones(new_eigval.n_rows, 1));
        }
        eigval -= eigval.min();
//        std::cout << "Eigval = " << std::endl << eigval << std::endl;
//        std::cout << "Degen = " << std::endl << degeneracy << std::endl;
//        std::cout << "S sqr = " << std::endl << s_squared_new_basis_vector << std::endl;
    }



private:
    const arma::dmat &js;
    unsigned int current_start;
    unsigned int current_end;

};

#endif //JULY_SUM_S_SQUARED_BLOCKED_H