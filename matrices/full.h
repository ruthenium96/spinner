#ifndef JULY_FULL_H
#define JULY_FULL_H

#include "base_matrices.h"

struct Full_Matrix_J : Matrices {
public:

    void add_to_s_squared(unsigned int i, unsigned int j, double v) override {
        s_squared(i, j) += v;
    };
    Add_To_Matrix s_squared_ = static_cast<Add_To_Matrix>(&Full_Matrix_J::add_to_s_squared);

    void add_to_hamiltonian(unsigned int i, unsigned int j, double v) override {
        hamiltonian(i, j) += v;
    }
    Add_To_Matrix hamiltonian_ = static_cast<Add_To_Matrix>(&Full_Matrix_J::add_to_hamiltonian);

    Full_Matrix_J(const Block_Lex_Map &blm,
                  const arma::dmat &js) : Matrices(blm), js(js) {
    }

    void construct_hamiltonian() {
        hamiltonian.resize(blm.tensor_size, blm.tensor_size);
#pragma omp parallel for
        for (unsigned long block = 0; block < blm.tensor_size; ++block) {
            for (int a = 0; a < blm.v_size; ++a) {
                for (int b = a + 1; b < blm.v_size; ++b) {
                    if (!std::isnan(js(a, b))) {
                        scalar_product_total(hamiltonian_, -2 * js(a, b), block, a, b);
                    }
                }
            }
        }
    }

    void construct_s_squared() {
        s_squared.resize(blm.tensor_size, blm.tensor_size);
        construct_block_of_s_squared(0, blm.tensor_size, s_squared_);
    }

    void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                            arma::vec & degeneracy) override{
        construct_hamiltonian();
        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, hamiltonian);
        eigval -= eigval.min();
        hamiltonian.reset();

        construct_s_squared();
        arma::mat s_squared_new_basis = eigvec.t() * s_squared * eigvec;
        s_squared_new_basis_vector = s_squared_new_basis.diag();
        degeneracy = arma::ones(size(eigval));
    }

private:
    const arma::dmat &js;
    arma::dmat hamiltonian;
    arma::dmat s_squared;
};

#endif //JULY_FULL_H
