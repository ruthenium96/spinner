//#ifndef JULY_SUM_BLOCKED_ONLY_MIDDLE_SPARSE_H
//#define JULY_SUM_BLOCKED_ONLY_MIDDLE_SPARSE_H
//
//#include "base_matrices.h"
//
//struct Only_J_Sparse : Matrices {
//public:
//    arma::sp_dmat middle_of_hamiltonian;
//    // it is a dense matrix!
//    // TODO: can I find eigenvalues of block_of_s_squared without building of matrix?
//    arma::dmat middle_of_s_squared;
//
//    void add_to_s_squared(unsigned int i, unsigned int j, double v) override {
//        middle_of_s_squared(i - start_of_middle, j - start_of_middle) += v;
//    };
//    Add_To_Matrix s_squared_ = static_cast<Add_To_Matrix>(&Only_J_Sparse::add_to_s_squared);
//
//    void add_to_hamiltonian(unsigned int i, unsigned int j, double v) override {
//        middle_of_hamiltonian(i - start_of_middle, j - start_of_middle) += v;
//    }
//    Add_To_Matrix hamiltonian_ = static_cast<Add_To_Matrix>(&Only_J_Sparse::add_to_hamiltonian);
//
//    Only_J_Sparse(const std::vector<int> &mults, const Indexes &blm,
//                  const arma::dmat &js) : Matrices(blm), js(js) {
//
//        unsigned long index = blm.block_boundaries.size() / 2 - 1;
//
//        start_of_middle = blm.block_boundaries[index];
//        end_of_middle = blm.block_boundaries[index + 1];
//
//        unsigned long size = end_of_middle - start_of_middle;
//        middle_of_hamiltonian.resize(size, size);
//        middle_of_s_squared.resize(size, size);
//    }
//
//    void construct_hamiltonian() {
////        #pragma omp parallel for
//        for (unsigned long block = start_of_middle; block < end_of_middle; ++block) {
//            for (int a = 0; a < blm.v_size; ++a) {
//                for (int b = a + 1; b < blm.v_size; ++b) {
//                    if (!std::isnan(js(a, b))) {
//                        scalar_product_total(hamiltonian_, -2 * js(a, b), block, a, b);
//                    }
//                }
//            }
//        }
////        std::cout  << "H:" << std::endl << middle_of_hamiltonian << std::endl;
//    }
//
//    void construct_s_squared() {
//        double diagonal_sum = 0;
//        for (double s : blm.mults) {
//            diagonal_sum += s * (s + 1);
//        }
////        #pragma omp parallel for
//        for (unsigned long block = start_of_middle; block < end_of_middle; ++block) {
//            add_to_s_squared(block, block, diagonal_sum);
//            for (int a = 0; a < blm.v_size; ++a) {
//                for (int b = a + 1; b < blm.v_size; ++b) {
//                    if (a != b) {
//                        scalar_product_total(s_squared_, 2, block, a, b);
//                    }
//                }
//            }
//        }
////        std::cout << "S^2:" << std::endl << middle_of_s_squared << std::endl;
//    }
//
//
//    void eigendecomposition() {
//        double g = 2.0;
//        // TODO: что-то вроде проверки на то, что матрица построена
//        arma::vec eigval;
//        arma::mat eigvec;
//        arma::eigs_sym(eigval, eigvec, middle_of_hamiltonian, 300, "sa");
//        eigval -= eigval.min();
//        arma::mat s_squared_new_basis = eigvec.t() * middle_of_s_squared * eigvec;
//        arma::vec s_squared_new_basis_vector = s_squared_new_basis.diag();
//        arma::vec degeneracities = sqrt(4 * s_squared_new_basis_vector + arma::ones(size(s_squared_new_basis_vector)));
////        std::cout << "Eigenvalues:" << std::endl << eigval << std::endl;
////        std::cout << "V^T * S^2 * V:" << std::endl << s_squared_new_basis_vector << std::endl;
////        std::cout << "sum_mults" << std::endl << degeneracities << std::endl;
//        arma::vec wep = s_squared_new_basis_vector % degeneracities;
////        std::cout << "wise-element product" << std::endl << wep << std::endl;
////        std::cout << "{" << std::endl;
//        for (int t = 1; t < 301; ++t) {
//            double mu_squared = g * g * sum(exp(- eigval / t) % wep);
//            double statsum = sum(exp(- eigval / t) % degeneracities);
////            std::cout << "T = " << t << ", chiT = " << mu_squared / statsum << std::endl;
////            std::cout << t << " : " << mu_squared / statsum << "," << std::endl;
//        }
////        std::cout << "}" << std::endl;
//    }
//
//private:
//    const arma::dmat &js;
//    unsigned int start_of_middle;
//    unsigned int end_of_middle;
//};
//
//#endif //JULY_SUM_BLOCKED_ONLY_MIDDLE_SPARSE_H
