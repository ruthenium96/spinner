#include "gtest/gtest.h"

#include "src/common/OneOrMany.h"

template<typename T>
T sum(T a, T b) {
    return a + b;
}

template<typename T>
T prod(T a, T b) {
    return a * b;
}

TEST(OneOrMany_tests, holdsOne) {
    OneOrMany<double> a = 12.0;
    EXPECT_TRUE(holdsOne(a));
    OneOrMany<double> b = std::vector<double>{12.0, 11.0, 10.0};
    EXPECT_FALSE(holdsOne(b));
}

TEST(OneOrMany_tests, holdsMany) {
    OneOrMany<double> a = 12.0;
    EXPECT_FALSE(holdsMany(a));
    OneOrMany<double> b = std::vector<double>{12.0, 11.0, 10.0};
    EXPECT_TRUE(holdsMany(b));
}


TEST(OneOrMany_tests, OneAndOneTransform) {
    OneOrMany<double> a = 2.0;
    OneOrMany<double> b = 3.0;

    OneOrMany<double> sum_ab = transform_one_or_many(std::function(sum<double>), a, b);
    OneOrMany<double> prod_ab = transform_one_or_many(std::function(prod<double>), a, b);
    ASSERT_TRUE(holdsOne(sum_ab));
    ASSERT_TRUE(holdsOne(prod_ab));
    EXPECT_EQ(getOneRef(sum_ab), 5.0);
    EXPECT_EQ(getOneRef(prod_ab), 6.0);
}

TEST(OneOrMany_tests, ManyAndManyTransform) {
    OneOrMany<double> a = std::vector<double>{2.0, 4.0, 6.0};
    OneOrMany<double> b = std::vector<double>{3.0, 5.0, 7.0};

    OneOrMany<double> sum_ab = transform_one_or_many(std::function(sum<double>), a, b);
    OneOrMany<double> prod_ab = transform_one_or_many(std::function(prod<double>), a, b);
    ASSERT_TRUE(holdsMany(sum_ab));
    ASSERT_TRUE(holdsMany(prod_ab));
    const auto& sum_ab_ref = getManyRef(sum_ab);
    const auto& prod_ab_ref = getManyRef(prod_ab);

    ASSERT_EQ(sum_ab_ref.size(), 3);
    ASSERT_EQ(prod_ab_ref.size(), 3);

    EXPECT_EQ(sum_ab_ref[0], 5.0);
    EXPECT_EQ(sum_ab_ref[1], 9.0);
    EXPECT_EQ(sum_ab_ref[2], 13.0);

    EXPECT_EQ(prod_ab_ref[0], 6.0);
    EXPECT_EQ(prod_ab_ref[1], 20.0);
    EXPECT_EQ(prod_ab_ref[2], 42.0);
}

TEST(OneOrMany_tests, ManyAndOneTransform) {
    OneOrMany<double> a = std::vector<double>{2.0, 4.0, 6.0};
    OneOrMany<double> b = 3.0;

    OneOrMany<double> sum_ab = transform_one_or_many(std::function(sum<double>), a, b);
    OneOrMany<double> prod_ab = transform_one_or_many(std::function(prod<double>), a, b);
    ASSERT_TRUE(holdsMany(sum_ab));
    ASSERT_TRUE(holdsMany(prod_ab));
    const auto& sum_ab_ref = getManyRef(sum_ab);
    const auto& prod_ab_ref = getManyRef(prod_ab);

    EXPECT_EQ(sum_ab_ref[0], 5.0);
    EXPECT_EQ(sum_ab_ref[1], 7.0);
    EXPECT_EQ(sum_ab_ref[2], 9.0);

    EXPECT_EQ(prod_ab_ref[0], 6.0);
    EXPECT_EQ(prod_ab_ref[1], 12.0);
    EXPECT_EQ(prod_ab_ref[2], 18.0);
}

TEST(OneOrMany_tests, OneAndManyTransform) {
    OneOrMany<double> a = 2.0;
    OneOrMany<double> b = std::vector<double>{3.0, 5.0, 7.0};

    OneOrMany<double> sum_ab = transform_one_or_many(std::function(sum<double>), a, b);
    OneOrMany<double> prod_ab = transform_one_or_many(std::function(prod<double>), a, b);
    ASSERT_TRUE(holdsMany(sum_ab));
    ASSERT_TRUE(holdsMany(prod_ab));
    const auto& sum_ab_ref = getManyRef(sum_ab);
    const auto& prod_ab_ref = getManyRef(prod_ab);

    ASSERT_EQ(sum_ab_ref.size(), 3);
    ASSERT_EQ(prod_ab_ref.size(), 3);

    EXPECT_EQ(sum_ab_ref[0], 5.0);
    EXPECT_EQ(sum_ab_ref[1], 7.0);
    EXPECT_EQ(sum_ab_ref[2], 9.0);

    EXPECT_EQ(prod_ab_ref[0], 6.0);
    EXPECT_EQ(prod_ab_ref[1], 10.0);
    EXPECT_EQ(prod_ab_ref[2], 14.0);
}