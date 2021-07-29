#include "indexes/blm.h"
#include "matrices/matrices.h"
#include <chrono>

// TODO: check if js_arma contains NaN

int main() {

    std::vector<int> mults = {2, 2, 2, 2};

    Indexes_P1 ind(mults);

    ind.construct_branching_diagram();
//
    arma::dmat js_arma(mults.size(), mults.size());
    // 0 - 2
    // 3 - 1
    double J = 10;
//    js_arma(0, 1) = J; js_arma(1, 0) = J;
    js_arma(0, 2) = J; js_arma(2, 0) = J;
    js_arma(1, 2) = J; js_arma(2, 1) = J;
    js_arma(1, 3) = J; js_arma(3, 1) = J;
    js_arma(3, 0) = J; js_arma(0, 3) = J;

//
    S2_Proj_Blocked_Matrices sm(ind, js_arma);
    arma::vec eigval;
    arma::vec s_squared_new_basis_vector;
    arma::vec degeneracy;
//
    sm.eigendecomposition(eigval, s_squared_new_basis_vector, degeneracy);

//    Full_Matrix_J fm(ind_typ, js_arma);
//    fm.eigendecomposition(eigval, s_squared_new_basis_vector, degeneracy);
//
//    std::cout << eigval << std::endl;



//    std::vector<int> mults = {2, 2, 2, 2, 2};
//    Block_Lex_Map_P2 blm(mults, 2);
//    blm.construct_branching_diagram();
//
//    arma::dmat js_arma(mults.size(), mults.size());
//    // 0 - 1
//    // 3 - 2
//    double J = 10;
//    js_arma(0, 1) = J; js_arma(1, 0) = J;
//    js_arma(1, 2) = J; js_arma(2, 1) = J;
//    js_arma(2, 3) = J; js_arma(3, 2) = J;
//    js_arma(3, 0) = J; js_arma(0, 3) = J;
//
//
//    Only_J_S2 jss(blm, js_arma);
//
//    arma::vec eigval;
//    arma::vec block_of_s_squared;
//    arma::vec degen;
//
//    jss.eigendecomposition(eigval, block_of_s_squared, degen);
//
//    std::cout << eigval << std::endl;

//    std::vector<int> mults = {4, 4, 4,
//                              4, 4, 4,
//                              4, 4, 4,};
//    Indexes blm(mults);
//    arma::dmat js_arma(mults.size(), mults.size());
//    // 0 - 1
//    // 3 - 2
//    js_arma(0, 1) = 50; js_arma(1, 0) = 50;
//    js_arma(1, 2) = 50; js_arma(2, 1) = 50;
//    js_arma(2, 3) = 50; js_arma(3, 2) = 50;
//    js_arma(3, 0) = 50; js_arma(0, 3) = 50;

//    // 0 - 1 - 2
//    // 5 - 4 - 3
//    js_arma(0, 1) = 5; js_arma(1, 0) = 5;
//    js_arma(1, 2) = 5; js_arma(2, 1) = 5;
//    js_arma(1, 4) = 5; js_arma(4, 1) = 5;
//    js_arma(2, 3) = 5; js_arma(3, 2) = 5;
//    js_arma(3, 4) = 5; js_arma(4, 3) = 5;
//    js_arma(4, 5) = 5; js_arma(5, 4) = 5;
//    js_arma(5, 0) = 5; js_arma(0, 5) = 5;

//    std::vector<int> mults = {4, 4, 4, 4,
//                              4, 4, 4, 4,};
//
////    std::vector<int> mults = {4, 4, 4, 4, 4,
////                              4, 4, 4, 4, 4,};
////
////
//    Indexes_P1 ind(mults);
//    ind.construct_branching_diagram();
//
//    arma::dmat js_arma(mults.size(), mults.size());
//
//    // 0 - 2 - 4 - 6
//    // 1 - 3 - 5 - 7
//    // vertical:
//    js_arma(0, 1) = 10; js_arma(1, 0) = 10;
//    js_arma(2, 3) = 10; js_arma(3, 2) = 10;
//    js_arma(4, 5) = 10; js_arma(5, 4) = 10;
//    js_arma(6, 7) = 10; js_arma(7, 6) = 10;
//
//    // horizontal:
//    js_arma(0, 2) = 10; js_arma(2, 0) = 10;
//    js_arma(2, 4) = 10; js_arma(4, 2) = 10;
//    js_arma(4, 6) = 10; js_arma(6, 4) = 10;
//    js_arma(1, 3) = 10; js_arma(3, 1) = 10;
//    js_arma(3, 5) = 10; js_arma(5, 3) = 10;
//    js_arma(5, 7) = 10; js_arma(7, 5) = 10;
////
////    // 0 - 2 - 4 - 6 - 8
////    // 1 - 3 - 5 - 7 - 9
////    // vertical:
////    js_arma(0, 1) = 10; js_arma(1, 0) = 10;
////    js_arma(2, 3) = 10; js_arma(3, 2) = 10;
////    js_arma(4, 5) = 10; js_arma(5, 4) = 10;
////    js_arma(6, 7) = 10; js_arma(7, 6) = 10;
////    js_arma(8, 9) = 10; js_arma(9, 8) = 10;
////
////    // horizontal:
////    js_arma(0, 2) = 10; js_arma(2, 0) = 10;
////    js_arma(2, 4) = 10; js_arma(4, 2) = 10;
////    js_arma(4, 6) = 10; js_arma(6, 4) = 10;
////    js_arma(6, 8) = 10; js_arma(8, 6) = 10;
////
////    js_arma(1, 3) = 10; js_arma(3, 1) = 10;
////    js_arma(3, 5) = 10; js_arma(5, 3) = 10;
////    js_arma(5, 7) = 10; js_arma(7, 5) = 10;
////    js_arma(7, 9) = 10; js_arma(9, 7) = 10;
//
//
////    std::vector<int> mults = {4, 4, 4,
////                              4, 4, 4,
////                              4, 4, 4,};
////
////    Indexes_P2 ind(mults, 4);
////
////    ind.construct_branching_diagram();
////
////    arma::dmat js_arma(mults.size(), mults.size());
////
////    // 0 - 2 - 4
////    // 7 - 8 - 6
////    // 5 - 3 - 1
////    // horizontal:
////    js_arma(0, 2) = 5; js_arma(2, 0) = 5;
////    js_arma(2, 4) = 5; js_arma(4, 2) = 5;
////    js_arma(7, 8) = 5; js_arma(8, 7) = 5;
////    js_arma(8, 6) = 5; js_arma(6, 8) = 5;
////    js_arma(5, 3) = 5; js_arma(3, 5) = 5;
////    js_arma(3, 1) = 5; js_arma(1, 3) = 5;
////
////    // vertical:
////    js_arma(0, 7) = 5; js_arma(7, 0) = 5;
////    js_arma(7, 5) = 5; js_arma(5, 7) = 5;
////    js_arma(2, 8) = 5; js_arma(8, 2) = 5;
////    js_arma(8, 3) = 5; js_arma(3, 8) = 5;
////    js_arma(4, 6) = 5; js_arma(6, 4) = 5;
////    js_arma(6, 1) = 5; js_arma(1, 6) = 5;
//
//    S2_Proj_Blocked_Matrices sm(ind, js_arma);
//    arma::vec eigval;
//    arma::vec s_squared_new_basis_vector;
//    arma::vec degeneracy;
//
//    auto start_dense = std::chrono::high_resolution_clock::now();
//    int cycles = 1;
//    for (int i = 0; i < cycles; ++i) {
//        sm.eigendecomposition(eigval, s_squared_new_basis_vector, degeneracy);
//    }
//    auto finish_dense = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed_dense_once = (finish_dense - start_dense) / cycles;
//    std::cout << "Elapsed time: " << elapsed_dense_once.count() << " s\n";

//    // 0 - 1 - 2
//    // 5 - 4 - 3
//    // 6 - 7 - 8
//    // horizontal:
//    js_arma(0, 1) = 5; js_arma(1, 0) = 5;
//    js_arma(1, 2) = 5; js_arma(2, 1) = 5;
//    js_arma(3, 4) = 5; js_arma(4, 3) = 5;
//    js_arma(4, 5) = 5; js_arma(5, 4) = 5;
//    js_arma(6, 7) = 5; js_arma(7, 6) = 5;
//    js_arma(7, 8) = 5; js_arma(8, 7) = 5;
//
//    // vertical:
//    js_arma(0, 5) = 5; js_arma(5, 0) = 5;
//    js_arma(5, 6) = 5; js_arma(6, 5) = 5;
//    js_arma(1, 4) = 5; js_arma(4, 1) = 5;
//    js_arma(4, 7) = 5; js_arma(7, 4) = 5;
//    js_arma(2, 3) = 5; js_arma(3, 2) = 5;
//    js_arma(3, 8) = 5; js_arma(8, 3) = 5;
//
//
//    Only_J ham(mults, blm, js_arma);
//    ham.construct_hamiltonian(); ham.construct_s_squared();
//    ham.eigendecomposition();

//    std::vector<int> mults = {2, 2,
//                              2, 2, 2, 2,
////                              2, 2
//    };
//    Indexes blm(mults);
//    arma::dmat js_arma(mults.size(), mults.size());
//    js_arma(0, 1) = 10; js_arma(1, 0) = 10;
//    js_arma(1, 2) = 10; js_arma(2, 1) = 10;
//    js_arma(2, 3) = 10; js_arma(3, 2) = 10;
//    js_arma(3, 4) = 10; js_arma(4, 3) = 10;
//    js_arma(4, 5) = 10; js_arma(5, 4) = 10;
////    js_arma(5, 6) = 10; js_arma(6, 5) = 10;
////    js_arma(6, 7) = 10; js_arma(7, 6) = 10;
//
//
////    Only_J_Sparse ham_sparse(mults, blm, js_arma);
////    ham_sparse.construct_hamiltonian(); ham_sparse.construct_s_squared();
////    ham_sparse.eigendecomposition();
//
//    Only_J ham_dense(mults, blm, js_arma);
//    ham_dense.construct_hamiltonian(); ham_dense.construct_s_squared();
//    ham_dense.eigendecomposition();


//    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
//                              2, 2,
//                              2, 2,
//                              //2
//                              };
//
//    arma::dmat js_arma(mults.size(), mults.size());
//    for (int i = 0; i < mults.size(); ++i) {
//        for (int j = i; j < mults.size(); ++j) {
//            js_arma(i, j) = NAN;
//        }
//    }
//    js_arma(0, 1) = 10; js_arma(1, 0) = 10;
//    js_arma(1, 2) = 10; js_arma(2, 1) = 10;
//    js_arma(2, 3) = 10; js_arma(3, 2) = 10;
//    js_arma(3, 4) = 10; js_arma(4, 3) = 10;
//    js_arma(4, 5) = 10; js_arma(5, 4) = 10;
//    js_arma(5, 6) = 10; js_arma(6, 5) = 10;
//    js_arma(6, 7) = 10; js_arma(7, 6) = 10;
//    js_arma(7, 8) = 10; js_arma(8, 7) = 10;
//    js_arma(8, 9) = 10; js_arma(9, 8) = 10;
//    js_arma(9, 10) = 10; js_arma(10, 9) = 10;
//    js_arma(10, 11) = 10; js_arma(11, 10) = 10;
//    js_arma(11, 12) = 10; js_arma(12, 11) = 10;
//    js_arma(12, 13) = 10; js_arma(13, 12) = 10;
////    js_arma(13, 14) = 140; js_arma(14, 13) = 140;
//
//    Indexes blm(mults);
//
//    auto start_dense = std::chrono::high_resolution_clock::now();
//    for (int i = 0; i < 10; ++i) {
//        Only_J ham(mults, blm, js_arma);
//        ham.construct_hamiltonian();
//        ham.construct_s_squared();
//        ham.eigendecomposition();
//    }
//    auto finish_dense = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed_dense = finish_dense - start_dense;
//    std::cout << "Elapsed time: " << elapsed_dense.count() << " s\n";
//
////    auto start_sparse = std::chrono::high_resolution_clock::now();
////    for (int i = 0; i < 10; ++i) {
////        Only_J_Sparse ham(mults, blm, js_arma);
////        ham.construct_hamiltonian();
////        ham.construct_s_squared();
////        ham.eigendecomposition();
////    }
////    auto finish_sparse = std::chrono::high_resolution_clock::now();
////    std::chrono::duration<double> elapsed_sparse = finish_sparse - start_sparse;
////    std::cout << "Elapsed time: " << elapsed_sparse.count() << " s\n";

//    std::vector<int> mults = {4, 4, 4, 4, 4, 4, 4, 4};
//
//    arma::dmat js_arma(mults.size(), mults.size());
//
//    js_arma(0, 2) = 10; js_arma(2, 0) = 10;
//    js_arma(2, 4) = 10; js_arma(4, 2) = 10;
//    js_arma(4, 6) = 10; js_arma(6, 4) = 10;
//    js_arma(6, 7) = 10; js_arma(7, 6) = 10;
//    js_arma(7, 5) = 10; js_arma(5, 7) = 10;
//    js_arma(5, 3) = 10; js_arma(3, 5) = 10;
//    js_arma(3, 1) = 10; js_arma(1, 3) = 10;
//
//
//    auto start = std::chrono::high_resolution_clock::now();
//    for (int i = 0; i < 0; ++i) {
//        Block_Lex_Map_P1 blm(mults);
//        blm.construct_branching_diagram();
//
//        Only_J_S2 jss(blm, js_arma);
//
//        arma::vec eigval;
//        arma::vec s_squared_new_basis_vector;
//        arma::vec degeneracy;
//        jss.eigendecomposition(eigval, s_squared_new_basis_vector, degeneracy);
//    }
//    auto finish = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = finish - start;
//    std::cout << "Elapsed time of P1: " << elapsed.count() << " s\n";
//
//    auto start_noopt = std::chrono::high_resolution_clock::now();
//    for (int i = 0; i < 1; ++i) {
//        Block_Lex_Map_P2 blm(mults, 4);
//        blm.construct_branching_diagram();
//
//        Only_J_S2 jss(blm, js_arma);
//
//        arma::vec eigval;
//        arma::vec s_squared_new_basis_vector;
//        arma::vec degeneracy;
//        jss.eigendecomposition(eigval, s_squared_new_basis_vector, degeneracy);
//    }
//    auto finish_noopt = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed_noopt = finish_noopt - start_noopt;
//    std::cout << "Elapsed time of P2: " << elapsed_noopt.count() << " s\n";
//
//    return 0;
}