#include "gtest/gtest.h"
#include <indexes/blm.h>

TEST(lex2block_correctness, 222) {
    std::vector<int> mults = {2, 2, 2};
    std::vector<int> correct_lex2block = {0, 1, 2, 4, 3, 5, 6, 7};
    Indexes_P1 blm(mults);
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        EXPECT_EQ(correct_lex2block[lex], blm.lex_to_block(lex));
    }
}

TEST(block2lex_correctness, 222) {
    std::vector<int> mults = {2, 2, 2};
    std::vector<int> correct_block2lex = {0, 1, 2, 4, 3, 5, 6, 7};
    Indexes_P1 blm(mults);
    for (int block = 0; block < blm.tensor_size; ++block) {
        EXPECT_EQ(correct_block2lex[block], blm.block_to_lex(block));
    }
}

TEST(lex2block_correctness, 2222) {
    std::vector<int> mults = {2, 2, 2, 2};
    std::vector<int> correct_lex2block = {0, 1, 2, 5, 3, 6, 7, 11, 4, 8, 9, 12, 10, 13, 14, 15};
    Indexes_P1 blm(mults);
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        EXPECT_EQ(correct_lex2block[lex], blm.lex_to_block(lex));
    }
}

TEST(block2lex_correctness, 2222) {
    std::vector<int> mults = {2, 2, 2, 2};
    std::vector<int> correct_block2lex = {0, 1, 2, 4, 8, 3, 5, 6, 9, 10, 12, 7, 11, 13, 14, 15};
    Indexes_P1 blm(mults);
    for (int block = 0; block < blm.tensor_size; ++block) {
        EXPECT_EQ(correct_block2lex[block], blm.block_to_lex(block));
    }
}

TEST(block_lex_reversibility, 2222) {
    std::vector<int> mults = {2, 2, 2, 2};
    Indexes_P1 blm(mults);
    for (int k = 0; k < blm.tensor_size; ++k) {
        EXPECT_EQ(k, blm.block_to_lex(blm.lex_to_block(k)));
        EXPECT_EQ(k, blm.lex_to_block(blm.block_to_lex(k)));
    }
}

TEST(block_lex_reversibility, 765432) {
    std::vector<int> mults = {7, 6, 5, 4, 3, 2};
    Indexes_P1 blm(mults);
    for (int k = 0; k < blm.tensor_size; ++k) {
        EXPECT_EQ(k, blm.block_to_lex(blm.lex_to_block(k)));
        EXPECT_EQ(k, blm.lex_to_block(blm.block_to_lex(k)));
    }
}

TEST(nzi_lex_reversibility, 2222) {
    std::vector<int> mults = {2, 2, 2, 2};
    Indexes_P1 blm(mults);
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        int calc_lex = 0;
        for (int i = 0; i < blm.v_size; ++i) {
            calc_lex += blm.nzi_to_lex(blm.lex_to_nzi(lex, i), i);
        }
        EXPECT_EQ(lex, calc_lex);
    }
}

TEST(nzi_lex_reversibility, 765432) {
    std::vector<int> mults = {7, 6, 5, 4, 3, 2};
    Indexes_P1 blm(mults);
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        int calc_lex = 0;
        for (int i = 0; i < blm.v_size; ++i) {
            calc_lex += blm.nzi_to_lex(blm.lex_to_nzi(lex, i), i);
        }
        EXPECT_EQ(lex, calc_lex);
    }
}

TEST(nzs_block_reversibility, 765432) {
    std::vector<int> mults = {7, 6, 5, 4, 3, 2};
    Indexes_P1 blm(mults);
    for (unsigned long block = 0; block < blm.tensor_size; ++block) {
        std::vector<int> nzs = blm.block_to_nzs(block);
        unsigned long calc_block = blm.nzs_to_block(nzs);
        EXPECT_EQ(block, calc_block);
    }
}


TEST(nzi_nzs_equivalence, 765432) {
    std::vector<int> mults = {7, 6, 5, 4, 3, 2};
    Indexes_P1 blm(mults);
    for (int block = 0; block < blm.tensor_size; ++block) {
        std::vector<int> nzs = blm.block_to_nzs(block);
        for (int i = 0; i < blm.v_size; ++i) {
            EXPECT_EQ(nzs[i], blm.block_to_nzi(block, i));
        }
    }
}

TEST(lexicographicity_of_lex, 765432) {
    std::vector<int> mults = {7, 6, 5, 4, 3, 2};
    Indexes_P1 blm(mults);
    std::vector<std::vector<int>> vectors(blm.tensor_size, std::vector<int>(blm.v_size, 0));
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        for (int i = 0; i < blm.v_size; ++i) {
            vectors[lex][i] = blm.lex_to_nzi(lex, i);
        }
        if (lex > 0) {
            ASSERT_TRUE(vectors[lex] > vectors[lex - 1]);
        }
    }
}

TEST(blocklexicographicity_of_block, 765432) {
    std::vector<int> mults = {7, 6, 5, 4, 3, 2};
    Indexes_P1 blm(mults);
    std::vector<std::vector<int>> vectors(blm.tensor_size, std::vector<int>(blm.v_size, 0));
    std::vector<int> sums(blm.tensor_size, 0);
    for (int block = 0; block < blm.tensor_size; ++block) {
        for (int i = 0; i < blm.v_size; ++i) {
            int proj = blm.block_to_nzi(block, i);
            vectors[block][i] = proj;
            sums[block] += proj;
        }
        if (block > 0) {
            bool lexicographicity_of_block = (sums[block] == sums[block - 1] && vectors[block] > vectors[block - 1]);
            bool start_of_new_block = (sums[block] == sums[block - 1] + 1);
            ASSERT_TRUE(lexicographicity_of_block || start_of_new_block);
        }
    }
}

TEST(lex_ladder_correctness, 222_zero) {
    std::vector<int> mults = {2, 2, 2};
    Indexes_P1 blm(mults);
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        for (int i = 0; i < blm.v_size; ++i) {
            int calc_lex = blm.lex_ladder(lex, i, 0);
            EXPECT_EQ(calc_lex, lex);
        }
    }
}

TEST(lex_ladder_correctness, 222_one) {
    std::vector<int> mults = {2, 2, 2};
    Indexes_P1 blm(mults);
    std::vector<int> nz(blm.tensor_size, 0);
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        // initialize current nz(lex)
        for (int i = 0; i < blm.v_size; ++i) {
            nz[i] = blm.lex_to_nzi(lex, i);
        }
        // for all centers we are trying ...
        for (int i = 0; i < blm.v_size; ++i) {
            // ... to decrease projection
            if (nz[i] > 0) {
                std::vector<int> nz_minus_one = nz;
                nz_minus_one[i] -= 1;
                int calc_lex = 0;
                // constructing calc_lex from nz
                for (int j = 0; j < blm.v_size; ++j) {
                    calc_lex += blm.nzi_to_lex(nz_minus_one[j], j);
                }
                EXPECT_EQ(calc_lex, blm.lex_ladder(lex, i, -1));
            }
            // ... to increase projection
            if (nz[i] < mults[i] - 1) {
                std::vector<int> nz_plus_one = nz;
                nz_plus_one[i] += 1;
                int calc_lex = 0;
                // constructing calc_lex from nz
                for (int j = 0; j < blm.v_size; ++j) {
                    calc_lex += blm.nzi_to_lex(nz_plus_one[j], j);
                }
                EXPECT_EQ(calc_lex, blm.lex_ladder(lex, i, 1));
            }
        }
    }
}

TEST(lex_ladder_correctness, 765432_one) {
    std::vector<int> mults = {7, 6, 5, 4, 3, 2};
    Indexes_P1 blm(mults);
    std::vector<int> nz(blm.tensor_size, 0);
    for (int lex = 0; lex < blm.tensor_size; ++lex) {
        // initialize current nz(lex)
        for (int i = 0; i < blm.v_size; ++i) {
            nz[i] = blm.lex_to_nzi(lex, i);
        }
        // for all centers we are trying ...
        for (int i = 0; i < blm.v_size; ++i) {
            // ... to decrease projection
            if (nz[i] > 0) {
                std::vector<int> nz_minus_one = nz;
                nz_minus_one[i] -= 1;
                int calc_lex = 0;
                // constructing calc_lex from nz
                for (int j = 0; j < blm.v_size; ++j) {
                    calc_lex += blm.nzi_to_lex(nz_minus_one[j], j);
                }
                EXPECT_EQ(calc_lex, blm.lex_ladder(lex, i, -1));
            }
            // ... to increase projection
            if (nz[i] < mults[i] - 1) {
                std::vector<int> nz_plus_one = nz;
                nz_plus_one[i] += 1;
                int calc_lex = 0;
                // constructing calc_lex from nz
                for (int j = 0; j < blm.v_size; ++j) {
                    calc_lex += blm.nzi_to_lex(nz_plus_one[j], j);
                }
                EXPECT_EQ(calc_lex, blm.lex_ladder(lex, i, 1));
            }
        }
    }
}

TEST(boundaries_correctness, _222) {
    std::vector<int> mults = {2, 2, 2};
    std::vector<int> correct_boundaries = {0, 1, 4, 7, 8};
    Indexes_P1 blm(mults);
    ASSERT_EQ(correct_boundaries.size(), blm.sym_sum_boundaries[0].size());
    for (int i = 0; i < correct_boundaries.size(); ++i) {
        EXPECT_EQ(correct_boundaries[i], blm.sym_sum_boundaries[0][i]);
    }
}

TEST(boundaries_correctness, _22) {
    std::vector<int> mults = {2, 2};
    std::vector<int> correct_boundaries = {0, 1, 3, 4};
    Indexes_P1 blm(mults);
    ASSERT_EQ(correct_boundaries.size(), blm.sym_sum_boundaries[0].size());
    for (int i = 0; i < correct_boundaries.size(); ++i) {
        EXPECT_EQ(correct_boundaries[i], blm.sym_sum_boundaries[0][i]);
    }
}