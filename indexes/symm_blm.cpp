//#include "symm_blm.h"
//
//void Block_Lex_Map_Symm::resize_CGs() {
//    boost::multi_array<double, 5>::extent_gen extents;
//    int max_t_spin_number = (int) (2 * max_t_spin + 1);
//    int max_n_spin_number = (int) (2 * max_n_spin + 1);
//    CGs.resize(extents[max_t_spin_number][max_n_spin_number][max_t_spin_number][2 * max_t_spin_number + 1][2 * max_n_spin_number + 1]);
//    std::fill_n(CGs.data(), CGs.num_elements(), NAN);
//}
//
//arma::dmat Block_Lex_Map_Symm::construct_transformation(int spin, int symmetry) const {
//    unsigned long br_start = spin_boundaries[spin].symm_boundaries[symmetry];
//    unsigned long br_end = spin_boundaries[spin].symm_boundaries[symmetry + 1];
//    unsigned long br_len = br_end - br_start;
//    unsigned long current_start = block_boundaries[spin];
//    unsigned long current_end = block_boundaries[spin + 1];
//    unsigned long bl_len = current_end - current_start;
//    arma::dmat trans(br_len, bl_len);
//    #pragma omp parallel for
//    for (unsigned long state = br_start; state < br_end; ++state) {
//        for (unsigned long block = current_start; block < current_end; ++block) {
//            trans(state - br_start, block - current_start) = total_coeff(state, block);
//        }
//    }
//
//    unsigned int full = 0;
//
//    for (int l = 0; l < trans.n_cols; ++l) {
//        for (int m = 0; m < trans.n_rows; ++m) {
//            if (trans(m, l) != 0) {
//                ++full;
//            }
//        }
//    }
//
//    std::cout << trans.n_rows << " " << trans.n_cols << " " << trans.n_rows * trans.n_cols * trans.n_cols << " " <<(double) full / (trans.n_cols * trans.n_rows) << std::endl;
//
//    return std::move(trans);
//}
//
//double Block_Lex_Map_Symm::hashed_clebsh_gordan(double l1, double l2, double l3, double m1,
//                                                double m2) const {
//    int il1 = (int) (2 * l1);
//    int il2 = (int) (2 * l2);
//    int il3 = (int) (2 * l3);
//    int im1 = (int) (2 * (m1 + max_t_spin));
//    int im2 = (int) (2 * (m2 + max_n_spin));
//    double& curr_value = CGs[il1][il2][il3][im1][im2];
//    if (std::isnan(curr_value)) {
//        curr_value = WignerSymbols::clebschGordan(l1, l2, l3,
//                                                  m1, m2,m1 + m2);
//    }
//    return curr_value;
//
//}
//
//Quantum_Numbers::Quantum_Numbers(int total_mult_, int symmetry_) : total_mult(total_mult_), symmetry(symmetry_) {}
