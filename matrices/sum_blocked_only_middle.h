#ifndef JULY_SUM_BLOCKED_ONLY_MIDDLE_H
#define JULY_SUM_BLOCKED_ONLY_MIDDLE_H

#include "base_matrices.h"

// build and diagonalize central block of Hamiltonian.
// then multiply formulas to multiplicity
struct Only_J : Matrices {
public:
    arma::dmat middle_of_hamiltonian;
    arma::dmat middle_of_s_squared;

    void add_to_s_squared(unsigned int i, unsigned int j, double v) override {
        middle_of_s_squared(i - start_of_middle, j - start_of_middle) += v;
    };
    Add_To_Matrix s_squared_ = static_cast<Add_To_Matrix>(&Only_J::add_to_s_squared);

    void add_to_hamiltonian(unsigned int i, unsigned int j, double v) override {
        middle_of_hamiltonian(i - start_of_middle, j - start_of_middle) += v;
    }
    Add_To_Matrix hamiltonian_ = static_cast<Add_To_Matrix>(&Only_J::add_to_hamiltonian);

    Only_J(const Block_Lex_Map &blm,
           const arma::dmat &js) : Matrices(blm), js(js) {

        unsigned long index = blm.block_boundaries.size() / 2 - 1;

        start_of_middle = blm.block_boundaries[index];
        end_of_middle = blm.block_boundaries[index + 1];

        unsigned long size = end_of_middle - start_of_middle;
        std::cout << size << std::endl;
        middle_of_hamiltonian.resize(size, size);
        middle_of_s_squared.resize(size, size);
    }

    void construct_hamiltonian() {
#pragma omp parallel for
        for (unsigned long block = start_of_middle; block < end_of_middle; ++block) {
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
        construct_block_of_s_squared(start_of_middle, end_of_middle, s_squared_);
    }

    void eigendecomposition(arma::vec & eigval, arma::vec & s_squared_new_basis_vector,
                            arma::vec & degeneracy) {
        // TODO: что-то вроде проверки на то, что матрица построена

        construct_hamiltonian();

        arma::mat eigvec;
        arma::eig_sym(eigval, eigvec, middle_of_hamiltonian);
        eigval -= eigval.min();
        middle_of_hamiltonian.reset();

        construct_s_squared();
        arma::mat s_squared_new_basis = eigvec.t() * middle_of_s_squared * eigvec;
        s_squared_new_basis_vector = s_squared_new_basis.diag();
        degeneracy = sqrt(4 * s_squared_new_basis_vector + arma::ones(size(s_squared_new_basis_vector)));
    }

private:
    const arma::dmat &js;
    unsigned int start_of_middle;
    unsigned int end_of_middle;
};

#endif //JULY_SUM_BLOCKED_ONLY_MIDDLE_H
