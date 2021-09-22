#include "gtest/gtest.h"
#include "entities/Space.h"
#include "components/Symmetrizer.h"
#include "components/TzSorter.h"
#include "components/NonAbelianSimplifier.h"
#include "common/Logger.h"

void compare_two_spaces(const Space& one, const Space& two) {
    EXPECT_EQ(one.blocks.size(), two.blocks.size())
    << "Sizes of spaces does not equal";
    for (size_t i = 0; i < one.blocks.size(); ++i) {
        int one_dimensionality = one.blocks[i].properties.dimensionality;
        int two_dimensionality = two.blocks[i].properties.dimensionality;
        int one_degeneracy = one.blocks[i].properties.degeneracy;
        int two_degeneracy = two.blocks[i].properties.degeneracy;

        EXPECT_EQ(one_dimensionality * one_degeneracy,
                  two_dimensionality * two_degeneracy)
                  << "First space :" << one.blocks[i].properties
                  << "Second space :" << two.blocks[i].properties;

        std::unordered_map<uint32_t, uint32_t> one_met;
        std::unordered_map<uint32_t, uint32_t> two_met;

        // counting:
        for (const auto& m : one.blocks[i].basis) {
            for (const auto p : m) {
                one_met[p.first] += 1;
            }
        }
        for (const auto& m : two.blocks[i].basis) {
            for (const auto p : m) {
                two_met[p.first] += 1;
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

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group(group::S3, {{1, 2, 0}, {0, 2, 1}});
    Symmetrizer symmetrizer(converter, group);

    space_full = symmetrizer.apply(space_full);
    space_simplified = symmetrizer.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333_doubleS3) {
    std::vector<int> mults = {3, 3, 3};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group(group::S3, {{1, 2, 0}, {0, 2, 1}});
    Symmetrizer symmetrizer(converter, group);

    space_full = symmetrizer.apply(space_full);
    space_full = symmetrizer.apply(space_full);
    space_simplified = symmetrizer.apply(space_simplified);
    space_simplified = symmetrizer.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333333_S3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group(group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    Symmetrizer symmetrizer(converter, group);

    space_full = symmetrizer.apply(space_full);
    space_simplified = symmetrizer.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333333_S3xS2) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group_S3(group::S3, {{1, 2, 0, 4, 5, 3}, {0, 2, 1, 3, 5, 4}});
    Symmetrizer symmetrizer_S3(converter, group_S3);

    Group group_S2(group::S2, {{3, 4, 5, 0, 1, 2}});
    Symmetrizer symmetrizer_S2(converter, group_S2);

    space_full = symmetrizer_S3.apply(space_full);
    space_full = symmetrizer_S2.apply(space_full);
    space_simplified = symmetrizer_S3.apply(space_simplified);
    space_simplified = symmetrizer_S2.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 222222222_after_second_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group_S3_first(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    Symmetrizer symmetrizer_S3_first(converter, group_S3_first);

    Group group_S3_second(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    Symmetrizer symmetrizer_S3_second(converter, group_S3_second);

    space_full = symmetrizer_S3_first.apply(space_full);
    space_full = symmetrizer_S3_second.apply(space_full);
    space_simplified = symmetrizer_S3_first.apply(space_simplified);
    space_simplified = symmetrizer_S3_second.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 222222222_after_first_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group_S3_first(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    Symmetrizer symmetrizer_S3_first(converter, group_S3_first);

    Group group_S3_second(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    Symmetrizer symmetrizer_S3_second(converter, group_S3_second);

    space_full = symmetrizer_S3_first.apply(space_full);
    space_full = symmetrizer_S3_second.apply(space_full);
    space_simplified = symmetrizer_S3_first.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    space_simplified = symmetrizer_S3_second.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 222222222_after_both_S3xS3) {
    std::vector<int> mults = {2, 2, 2, 2, 2, 2, 2, 2, 2};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group_S3_first(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    Symmetrizer symmetrizer_S3_first(converter, group_S3_first);

    Group group_S3_second(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    Symmetrizer symmetrizer_S3_second(converter, group_S3_second);

    NonAbelianSimplifier nonAbelianSimplifier;

    space_full = symmetrizer_S3_first.apply(space_full);
    space_full = symmetrizer_S3_second.apply(space_full);

    space_simplified = symmetrizer_S3_first.apply(space_simplified);
    space_simplified = nonAbelianSimplifier.apply(space_simplified);
    space_simplified = symmetrizer_S3_second.apply(space_simplified);
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333333333_after_second_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group_S3_first(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    Symmetrizer symmetrizer_S3_first(converter, group_S3_first);

    Group group_S3_second(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    Symmetrizer symmetrizer_S3_second(converter, group_S3_second);

    space_full = symmetrizer_S3_first.apply(space_full);
    space_full = symmetrizer_S3_second.apply(space_full);
    space_simplified = symmetrizer_S3_first.apply(space_simplified);
    space_simplified = symmetrizer_S3_second.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333333333_after_first_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3,};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group_S3_first(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    Symmetrizer symmetrizer_S3_first(converter, group_S3_first);

    Group group_S3_second(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    Symmetrizer symmetrizer_S3_second(converter, group_S3_second);

    space_full = symmetrizer_S3_first.apply(space_full);
    space_full = symmetrizer_S3_second.apply(space_full);
    space_simplified = symmetrizer_S3_first.apply(space_simplified);

    NonAbelianSimplifier nonAbelianSimplifier;
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    space_simplified = symmetrizer_S3_second.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}

TEST(nonAbelianSimplifier, 333333333_after_both_S3xS3) {
    std::vector<int> mults = {3, 3, 3, 3, 3, 3, 3, 3, 3,};

    spaces::LexicographicIndexConverter converter(mults);

    Space space_full(converter.total_space_size);
    Space space_simplified(converter.total_space_size);

    Group group_S3_first(group::S3, {{1, 2, 0, 4, 5, 3, 7, 8, 6}, {0, 2, 1, 3, 5, 4, 6, 8, 7}});
    Symmetrizer symmetrizer_S3_first(converter, group_S3_first);

    Group group_S3_second(group::S3, {{3, 4, 5, 6, 7, 8, 0, 1, 2}, {0, 1, 2, 6, 7, 8, 3, 4, 5}});
    Symmetrizer symmetrizer_S3_second(converter, group_S3_second);

    NonAbelianSimplifier nonAbelianSimplifier;

    space_full = symmetrizer_S3_first.apply(space_full);
    space_full = symmetrizer_S3_second.apply(space_full);

    space_simplified = symmetrizer_S3_first.apply(space_simplified);
    space_simplified = nonAbelianSimplifier.apply(space_simplified);
    space_simplified = symmetrizer_S3_second.apply(space_simplified);
    space_simplified = nonAbelianSimplifier.apply(space_simplified);

    compare_two_spaces(space_full, space_simplified);
}