#include "gtest/gtest.h"
#include "common/Runner.h"
#include "common/Logger.h"

#include <cmath>

void compare_two_spaces(const Space& one, const Space& two) {
    EXPECT_EQ(one.blocks.size(), two.blocks.size())
    << "Sizes of spaces does not equal";
    for (size_t i = 0; i < one.blocks.size(); ++i) {
        uint32_t one_dimensionality = one.blocks[i].properties.dimensionality;
        uint32_t two_dimensionality = two.blocks[i].properties.dimensionality;
        uint32_t one_degeneracy = one.blocks[i].properties.degeneracy;
        uint32_t two_degeneracy = two.blocks[i].properties.degeneracy;

        EXPECT_EQ(one_dimensionality * one_degeneracy,
                  two_dimensionality * two_degeneracy)
                  << "First space :" << one.blocks[i].properties
                  << "Second space :" << two.blocks[i].properties;

        std::unordered_map<uint32_t, uint32_t> one_met;
        std::unordered_map<uint32_t, uint32_t> two_met;

        // counting:
        for (uint32_t index_of_vector = 0; index_of_vector < one.blocks[i].decomposition.size(); ++index_of_vector) {
            auto iterator = one.blocks[i].decomposition.GetNewIterator(index_of_vector);
            while (iterator->hasNext()) {
                auto item = iterator->getNext();
                one_met[item.index] += 1;
            }
        }
        for (uint32_t index_of_vector = 0; index_of_vector < two.blocks[i].decomposition.size(); ++index_of_vector) {
            auto iterator = two.blocks[i].decomposition.GetNewIterator(index_of_vector);
            while (iterator->hasNext()) {
                auto item = iterator->getNext();
                two_met[item.index] += 1;
            }
        }

        // ceiling to closest n * dimensionality, because some projectors has some coefficients equal to 0
        for (auto& p : one_met) {
            p.second = (int) ceil(((double) p.second) / one_dimensionality) * one_dimensionality;
        }
        for (auto& p : two_met) {
            p.second = (int) ceil(((double) p.second) / two_dimensionality) * two_dimensionality;
        }

        // comparing spaces
        for (const auto& p : one_met) {
            EXPECT_EQ(p.second * one_degeneracy,two_met[p.first] * two_degeneracy)
                      << "Lex index [" << p.first << "] is met different times in spaces, representation: "
                      << one.blocks[i].properties.get_representation_name();
        }
    }
}

TEST(nonAbelianSimplifier, 333_S3) {
    std::vector<int> mults = {3, 3, 3};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner_simplified.NonAbelianSimplify();

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 333_doubleS3) {

    std::vector<int> mults = {3, 3, 3};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner_full.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner_simplified.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner_simplified.NonAbelianSimplify();

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 333_doubleS3_tricky) {

    std::vector<int> mults = {3, 3, 3};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner_full.Symmetrize(group::S3, {{2, 0, 1}, {1, 0, 2}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0}, {0, 2, 1}});
    runner_simplified.Symmetrize(group::S3, {{2, 0, 1}, {1, 0, 2}});
    runner_simplified.NonAbelianSimplify();

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 333333_S3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    runner_simplified.NonAbelianSimplify();

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 333333_S3xS2) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    runner_full.Symmetrize(group::S2, {{3, 4, 5, 0, 1, 2}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    runner_simplified.Symmetrize(group::S2, {{3, 4, 5, 0, 1, 2}});
    runner_simplified.NonAbelianSimplify();

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 222222222_after_second_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_full.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_simplified.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    EXPECT_THROW(runner_simplified.NonAbelianSimplify(), std::invalid_argument);

//    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 222222222_after_first_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_full.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_simplified.NonAbelianSimplify();
    runner_simplified.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 222222222_after_both_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_full.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_simplified.NonAbelianSimplify();
    runner_simplified.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    runner_simplified.NonAbelianSimplify();

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 333333333_after_second_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_full.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_simplified.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    EXPECT_THROW(runner_simplified.NonAbelianSimplify(), std::invalid_argument);

//    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 333333333_after_first_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3,};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_full.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_simplified.NonAbelianSimplify();
    runner_simplified.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}

TEST(nonAbelianSimplifier, 333333333_after_both_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3,};

    runner::Runner runner_full(mults);
    runner::Runner runner_simplified(mults);

    runner_full.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_full.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});

    runner_simplified.Symmetrize(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    runner_simplified.NonAbelianSimplify();
    runner_simplified.Symmetrize(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    runner_simplified.NonAbelianSimplify();

    compare_two_spaces(runner_full.getSpace(), runner_simplified.getSpace());
}