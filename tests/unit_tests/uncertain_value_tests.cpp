#include <cstddef>
#include <stdexcept>
#include "gtest/gtest.h"
#include "magic_enum.hpp"
#include "src/common/UncertainValue.h"

TEST(UncertainValueTest, DefaultConstructor) {
    common::UncertainValue a;
    EXPECT_DOUBLE_EQ(a.mean(), 0.0);
    for (size_t i = 0; i < magic_enum::enum_count<common::UncertaintySources>(); ++i) {
        EXPECT_DOUBLE_EQ(a.stdevs()[i], 0.0);
    }
    EXPECT_DOUBLE_EQ(a.stdev_total(), 0.0);
}

TEST(UncertainValueTest, SingleValueConstructor) {
    common::UncertainValue a(5.0);
    EXPECT_DOUBLE_EQ(a.mean(), 5.0);
    for (size_t i = 0; i < magic_enum::enum_count<common::UncertaintySources>(); ++i) {
        EXPECT_DOUBLE_EQ(a.stdevs()[i], 0.0);
    }
    EXPECT_DOUBLE_EQ(a.stdev_total(), 0.0);
}

TEST(UncertainValueTest, SourceSpecificConstructor) {
    common::UncertainValue a(10.0, 0.5, common::UncertaintySources::FTLM);
    EXPECT_DOUBLE_EQ(a.mean(), 10.0);
    EXPECT_DOUBLE_EQ(a.stdevs()[common::UncertaintySources::FTLM], 0.5);
    EXPECT_DOUBLE_EQ(a.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, SourceSpecificConstructorThrowNegativeStdev) {
    EXPECT_THROW(common::UncertainValue(10.0, -0.0005, common::UncertaintySources::FTLM), std::invalid_argument);
}

TEST(UncertainValueTest, ArrayConstructor) {
    common::UncertainValue a(3.0, {0.1, 0.2});
    EXPECT_DOUBLE_EQ(a.mean(), 3.0);
    EXPECT_DOUBLE_EQ(a.stdevs()[common::UncertaintySources::FTLM], 0.1);
    EXPECT_DOUBLE_EQ(a.stdevs()[common::UncertaintySources::EXPERIMENT], 0.2);
    EXPECT_DOUBLE_EQ(a.stdev_total(), std::sqrt(0.01 + 0.04));
}

TEST(UncertainValueTest, ArrayConstructorThrowNegativeStdev) {
    EXPECT_THROW(common::UncertainValue(3.0, {0.1, -0.00002}), std::invalid_argument);
}

TEST(UncertainValueTest, AdditionBasicFTLM) {
    common::UncertainValue a(10.0, 0.3, common::UncertaintySources::FTLM);
    common::UncertainValue b(15.0, 0.4, common::UncertaintySources::FTLM);
    common::UncertainValue c = a + b;
    EXPECT_DOUBLE_EQ(c.mean(), 25.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 0.7); // 0.3 + 0.4
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
    common::UncertainValue d = a;
    d += b;
    EXPECT_DOUBLE_EQ(d.mean(), 25.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::FTLM], 0.7); // 0.3 + 0.4
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, AdditionBasicEXPERIMENT) {
    common::UncertainValue a(10.0, 0.3, common::UncertaintySources::EXPERIMENT);
    common::UncertainValue b(15.0, 0.4, common::UncertaintySources::EXPERIMENT);
    common::UncertainValue c = a + b;
    EXPECT_DOUBLE_EQ(c.mean(), 25.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], std::sqrt(0.09 + 0.16));
    common::UncertainValue d = a;
    d += b;
    EXPECT_DOUBLE_EQ(d.mean(), 25.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::EXPERIMENT], std::sqrt(0.09 + 0.16));
}

TEST(UncertainValueTest, AdditionDifferentSources) {
    common::UncertainValue a(10.0, 0.3, common::UncertaintySources::FTLM);
    common::UncertainValue b(15.0, 0.4, common::UncertaintySources::EXPERIMENT);
    common::UncertainValue c = a + b;
    EXPECT_DOUBLE_EQ(c.mean(), 25.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], 0.4);
    EXPECT_DOUBLE_EQ(c.stdev_total(), std::sqrt(0.09 + 0.16));
    common::UncertainValue d = a;
    d += b;
    EXPECT_DOUBLE_EQ(d.mean(), 25.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::EXPERIMENT], 0.4);
    EXPECT_DOUBLE_EQ(d.stdev_total(), std::sqrt(0.09 + 0.16));
}

TEST(UncertainValueTest, AdditionInverse) {
    common::UncertainValue a(10.0, {0.3, 0.4});
    common::UncertainValue b = a + (-a);
    EXPECT_DOUBLE_EQ(b.mean(), 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::UncertaintySources::EXPERIMENT], std::sqrt(0.32));
    EXPECT_DOUBLE_EQ(b.stdev_total(), b.stdevs()[common::UncertaintySources::EXPERIMENT]);
}

TEST(UncertainValueTest, AdditionCommutativity) {
    common::UncertainValue a(10.0, {0.3, 0.4});
    common::UncertainValue b(-15.0, {0.4, 0.5});
    common::UncertainValue c = a + b;
    common::UncertainValue d = b + a;
    EXPECT_DOUBLE_EQ(c.mean(), d.mean());
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], d.stdevs()[common::UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], d.stdevs()[common::UncertaintySources::EXPERIMENT]);
    EXPECT_DOUBLE_EQ(c.stdev_total(), d.stdev_total());
}

TEST(UncertainValueTest, AdditionAssociativity) {
    common::UncertainValue a(10.0, {0.3, 0.4});
    common::UncertainValue b(-15.0, {0.4, 0.5});
    common::UncertainValue c(1.0, {0.5, 0.6});
    common::UncertainValue d = (a + b) + c;
    common::UncertainValue e = a + (b + c);
    EXPECT_DOUBLE_EQ(d.mean(), e.mean());
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::FTLM], e.stdevs()[common::UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::EXPERIMENT], e.stdevs()[common::UncertaintySources::EXPERIMENT]);
    EXPECT_DOUBLE_EQ(d.stdev_total(), e.stdev_total());
}

TEST(UncertainValueTest, SubtractionFTLM) {
    common::UncertainValue a(20.0, 0.5, common::UncertaintySources::FTLM);
    common::UncertainValue b(15.0, 0.2, common::UncertaintySources::FTLM);
    common::UncertainValue c = a - b;
    EXPECT_DOUBLE_EQ(c.mean(), 5.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
    common::UncertainValue d = a;
    d -= b;
    EXPECT_DOUBLE_EQ(d.mean(), 5.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, SubtractionEXPERIMENT) {
    common::UncertainValue a(20.0, 0.5, common::UncertaintySources::EXPERIMENT);
    common::UncertainValue b(15.0, 0.2, common::UncertaintySources::EXPERIMENT);
    common::UncertainValue c = a - b;
    EXPECT_DOUBLE_EQ(c.mean(), 5.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], std::sqrt(0.25+ 0.04));
    common::UncertainValue d = a;
    d -= b;
    EXPECT_DOUBLE_EQ(d.mean(), 5.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::EXPERIMENT], std::sqrt(0.25+ 0.04));
}

TEST(UncertainValueTest, SubtractionAdditionOfUnaryNegatedFTLM) {
    common::UncertainValue a(20.0, 0.5, common::UncertaintySources::FTLM);
    common::UncertainValue b(15.0, 0.25, common::UncertaintySources::FTLM);
    common::UncertainValue c = a - b;
    common::UncertainValue d = a + (-b);
    EXPECT_DOUBLE_EQ(c.mean(), d.mean());
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], d.stdevs()[common::UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], d.stdevs()[common::UncertaintySources::EXPERIMENT]);
    common::UncertainValue e = a + (-1.0 * b);
    EXPECT_DOUBLE_EQ(c.mean(), e.mean());
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], e.stdevs()[common::UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], e.stdevs()[common::UncertaintySources::EXPERIMENT]);
    common::UncertainValue f = b - a;
    common::UncertainValue g = b + (-a);
    EXPECT_DOUBLE_EQ(f.mean(), g.mean());
    EXPECT_DOUBLE_EQ(f.stdevs()[common::UncertaintySources::FTLM], g.stdevs()[common::UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(f.stdevs()[common::UncertaintySources::EXPERIMENT], g.stdevs()[common::UncertaintySources::EXPERIMENT]);
    common::UncertainValue h = b + (-1.0 * a);
    EXPECT_DOUBLE_EQ(f.mean(), h.mean());
    EXPECT_DOUBLE_EQ(f.stdevs()[common::UncertaintySources::FTLM], h.stdevs()[common::UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(f.stdevs()[common::UncertaintySources::EXPERIMENT], h.stdevs()[common::UncertaintySources::EXPERIMENT]);
}

TEST(UncertainValueTest, SubtractionCorrelationSignHandlingFTLM) {
    common::UncertainValue a(20.0, 0.3, common::UncertaintySources::FTLM);
    common::UncertainValue b(15.0, 0.5, common::UncertaintySources::FTLM);
    common::UncertainValue c(10.0, 0.06, common::UncertaintySources::FTLM);
    common::UncertainValue d = (a - b) + c;
    EXPECT_DOUBLE_EQ(d.mean(), 15.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::FTLM], 0.14);
    EXPECT_DOUBLE_EQ(d.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
    common::UncertainValue e = c + (a - b);
    EXPECT_DOUBLE_EQ(e.mean(), 15.0);
    EXPECT_DOUBLE_EQ(e.stdevs()[common::UncertaintySources::FTLM], 0.14);
    EXPECT_DOUBLE_EQ(e.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
    common::UncertainValue f = c - (b - a);
    EXPECT_DOUBLE_EQ(f.mean(), 15.0);
    EXPECT_DOUBLE_EQ(f.stdevs()[common::UncertaintySources::FTLM], 0.14);
    EXPECT_DOUBLE_EQ(f.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);

    common::UncertainValue g = (a - b) + c * 5;
    EXPECT_DOUBLE_EQ(g.mean(), 55.0);
    EXPECT_NEAR(g.stdevs()[common::UncertaintySources::FTLM], 0.10, 1e-6);
    EXPECT_DOUBLE_EQ(g.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, MultiplicationBasicFTLM) {
    common::UncertainValue a(3.0, 0.3, common::UncertaintySources::FTLM);
    common::UncertainValue b(4.0, 0.4, common::UncertaintySources::FTLM);
    common::UncertainValue c = a * b;
    EXPECT_DOUBLE_EQ(c.mean(), 12.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 2.4);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, MultiplicationBasicEXPERIMENT) {
    common::UncertainValue a(3.0, 0.3, common::UncertaintySources::EXPERIMENT);
    common::UncertainValue b(4.0, 0.4, common::UncertaintySources::EXPERIMENT);
    common::UncertainValue c = a * b;
    EXPECT_DOUBLE_EQ(c.mean(), 12.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::EXPERIMENT], 12*std::sqrt(0.02));
}

TEST(UncertainValueTest, MultiplicationNegativeFTLM) {
    common::UncertainValue a(5.0, 0.5, common::UncertaintySources::FTLM);
    common::UncertainValue b(-2.0, 0.2, common::UncertaintySources::FTLM);
    common::UncertainValue c = a * b;
    EXPECT_DOUBLE_EQ(c.mean(), -10.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[common::UncertaintySources::FTLM], 2.0);
}

TEST(UncertainValueTest, MultiplicationInverse) {
    common::UncertainValue a(10.0, {0.3, 0.4});
    common::UncertainValue b = a / a;
    EXPECT_DOUBLE_EQ(b.mean(), 1.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::UncertaintySources::EXPERIMENT], std::sqrt(0.0032));
    EXPECT_DOUBLE_EQ(b.stdev_total(), b.stdevs()[common::UncertaintySources::EXPERIMENT]);
}

TEST(UncertainValueTest, MultiplyByZero) {
    common::UncertainValue a(10.0, {0.3, 0.4});
    common::UncertainValue zero(0.0);
    common::UncertainValue b = a * zero;
    EXPECT_DOUBLE_EQ(b.mean(), 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, ScalarMultiplication) {
    common::UncertainValue a(10.0, {0.3, 0.4});
    common::UncertainValue b = 2.0 * a;
    EXPECT_DOUBLE_EQ(b.mean(), 20.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::FTLM], 0.6);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::EXPERIMENT], 0.8);
}

TEST(UncertainValueTest, SqrtBasic) {
    common::UncertainValue a(16.0, 1.6, common::UncertaintySources::FTLM);
    common::UncertainValue b = common::UncertainValue::sqrt(a);
    EXPECT_DOUBLE_EQ(b.mean(), 4.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[common::UncertaintySources::FTLM], 0.2);
}

TEST(UncertainValueTest, SqrtSquareInversibilityFTLM) {
    common::UncertainValue a(16.0, 0.25, common::UncertaintySources::FTLM);
    common::UncertainValue b = common::UncertainValue::sqrt(a);
    common::UncertainValue c = b * b;
    EXPECT_DOUBLE_EQ(a.mean(), c.mean());
    EXPECT_DOUBLE_EQ(a.stdevs()[common::UncertaintySources::FTLM], c.stdevs()[common::UncertaintySources::FTLM]);
}

TEST(UncertainValueTest, SquareSqrtInversibilityFTLM) {
    common::UncertainValue a(16.0, 0.25, common::UncertaintySources::FTLM);
    common::UncertainValue b = a * a;
    common::UncertainValue c = common::UncertainValue::sqrt(b);
    EXPECT_DOUBLE_EQ(a.mean(), c.mean());
    EXPECT_DOUBLE_EQ(a.stdevs()[common::UncertaintySources::FTLM], c.stdevs()[common::UncertaintySources::FTLM]);
}

TEST(UncertainValueTest, InverseIdentity) {
    common::UncertainValue a(5.0, {0.1, 0.2});
    common::UncertainValue b = common::UncertainValue::inv(common::UncertainValue::inv(a));
    EXPECT_NEAR(a.mean(), b.mean(), 1e-10);
    EXPECT_NEAR(a.stdevs()[common::FTLM], b.stdevs()[common::FTLM], 1e-10);
    EXPECT_NEAR(a.stdevs()[common::EXPERIMENT], b.stdevs()[common::EXPERIMENT], 1e-10);
}

TEST(UncertainValueTest, MultiplicationSignHandling) {
    common::UncertainValue pos(5.0, 0.5, common::UncertaintySources::FTLM);
    common::UncertainValue neg(-3.0, 0.3, common::UncertaintySources::FTLM);
    
    common::UncertainValue pn = pos * neg;
    EXPECT_DOUBLE_EQ(pn.mean(), -15.0);
    EXPECT_NEAR(pn.stdevs()[common::UncertaintySources::FTLM], 
                3.0, 
                1e-10);

    common::UncertainValue nn = neg * neg;
    EXPECT_DOUBLE_EQ(nn.mean(), 9.0);
    EXPECT_NEAR(nn.stdevs()[common::UncertaintySources::FTLM], 
                1.8, 
                1e-10);
}

TEST(UncertainValueTest, DivisionFTLM) {
    common::UncertainValue a(10.0, 0.5, common::UncertaintySources::FTLM);
    common::UncertainValue b(-2.0, 0.1, common::UncertaintySources::FTLM);

    common::UncertainValue c = a / b;
    EXPECT_DOUBLE_EQ(c.mean(), -5.0);

    EXPECT_NEAR(c.stdevs()[common::UncertaintySources::FTLM], 0.0, 1e-10);
}