#include <cstddef>
#include <stdexcept>
#include "gtest/gtest.h"
#include "magic_enum.hpp"
#include "src/common/UncertainValue.h"

using namespace spinner::common;

TEST(UncertainValueTest, DefaultConstructor) {
    UncertainValue a;
    EXPECT_DOUBLE_EQ(a.mean(), 0.0);
    for (size_t i = 0; i < magic_enum::enum_count<UncertaintySources>(); ++i) {
        EXPECT_DOUBLE_EQ(a.stdevs()[i], 0.0);
    }
    EXPECT_DOUBLE_EQ(a.stdev_total(), 0.0);
}

TEST(UncertainValueTest, SingleValueConstructor) {
    UncertainValue a(5.0);
    EXPECT_DOUBLE_EQ(a.mean(), 5.0);
    for (size_t i = 0; i < magic_enum::enum_count<UncertaintySources>(); ++i) {
        EXPECT_DOUBLE_EQ(a.stdevs()[i], 0.0);
    }
    EXPECT_DOUBLE_EQ(a.stdev_total(), 0.0);
}

TEST(UncertainValueTest, SourceSpecificConstructor) {
    UncertainValue a(10.0, 0.5, UncertaintySources::FTLM);
    EXPECT_DOUBLE_EQ(a.mean(), 10.0);
    EXPECT_DOUBLE_EQ(a.stdevs()[UncertaintySources::FTLM], 0.5);
    EXPECT_DOUBLE_EQ(a.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, SourceSpecificConstructorThrowNegativeStdev) {
    EXPECT_THROW(UncertainValue(10.0, -0.0005, UncertaintySources::FTLM), std::invalid_argument);
}

TEST(UncertainValueTest, ArrayConstructor) {
    UncertainValue a(3.0, {0.1, 0.2});
    EXPECT_DOUBLE_EQ(a.mean(), 3.0);
    EXPECT_DOUBLE_EQ(a.stdevs()[UncertaintySources::FTLM], 0.1);
    EXPECT_DOUBLE_EQ(a.stdevs()[UncertaintySources::EXPERIMENT], 0.2);
    EXPECT_DOUBLE_EQ(a.stdev_total(), std::sqrt(0.01 + 0.04));
}

TEST(UncertainValueTest, ArrayConstructorThrowNegativeStdev) {
    EXPECT_THROW(UncertainValue(3.0, {0.1, -0.00002}), std::invalid_argument);
}

TEST(UncertainValueTest, AdditionBasicFTLM) {
    UncertainValue a(10.0, 0.3, UncertaintySources::FTLM);
    UncertainValue b(15.0, 0.4, UncertaintySources::FTLM);
    UncertainValue c = a + b;
    EXPECT_DOUBLE_EQ(c.mean(), 25.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 0.7); // 0.3 + 0.4
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
    UncertainValue d = a;
    d += b;
    EXPECT_DOUBLE_EQ(d.mean(), 25.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::FTLM], 0.7); // 0.3 + 0.4
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, AdditionBasicEXPERIMENT) {
    UncertainValue a(10.0, 0.3, UncertaintySources::EXPERIMENT);
    UncertainValue b(15.0, 0.4, UncertaintySources::EXPERIMENT);
    UncertainValue c = a + b;
    EXPECT_DOUBLE_EQ(c.mean(), 25.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], std::sqrt(0.09 + 0.16));
    UncertainValue d = a;
    d += b;
    EXPECT_DOUBLE_EQ(d.mean(), 25.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::EXPERIMENT], std::sqrt(0.09 + 0.16));
}

TEST(UncertainValueTest, AdditionDifferentSources) {
    UncertainValue a(10.0, 0.3, UncertaintySources::FTLM);
    UncertainValue b(15.0, 0.4, UncertaintySources::EXPERIMENT);
    UncertainValue c = a + b;
    EXPECT_DOUBLE_EQ(c.mean(), 25.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], 0.4);
    EXPECT_DOUBLE_EQ(c.stdev_total(), std::sqrt(0.09 + 0.16));
    UncertainValue d = a;
    d += b;
    EXPECT_DOUBLE_EQ(d.mean(), 25.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::EXPERIMENT], 0.4);
    EXPECT_DOUBLE_EQ(d.stdev_total(), std::sqrt(0.09 + 0.16));
}

TEST(UncertainValueTest, AdditionInverse) {
    UncertainValue a(10.0, {0.3, 0.4});
    UncertainValue b = a + (-a);
    EXPECT_DOUBLE_EQ(b.mean(), 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[UncertaintySources::EXPERIMENT], std::sqrt(0.32));
    EXPECT_DOUBLE_EQ(b.stdev_total(), b.stdevs()[UncertaintySources::EXPERIMENT]);
}

TEST(UncertainValueTest, AdditionCommutativity) {
    UncertainValue a(10.0, {0.3, 0.4});
    UncertainValue b(-15.0, {0.4, 0.5});
    UncertainValue c = a + b;
    UncertainValue d = b + a;
    EXPECT_DOUBLE_EQ(c.mean(), d.mean());
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], d.stdevs()[UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], d.stdevs()[UncertaintySources::EXPERIMENT]);
    EXPECT_DOUBLE_EQ(c.stdev_total(), d.stdev_total());
}

TEST(UncertainValueTest, AdditionAssociativity) {
    UncertainValue a(10.0, {0.3, 0.4});
    UncertainValue b(-15.0, {0.4, 0.5});
    UncertainValue c(1.0, {0.5, 0.6});
    UncertainValue d = (a + b) + c;
    UncertainValue e = a + (b + c);
    EXPECT_DOUBLE_EQ(d.mean(), e.mean());
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::FTLM], e.stdevs()[UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::EXPERIMENT], e.stdevs()[UncertaintySources::EXPERIMENT]);
    EXPECT_DOUBLE_EQ(d.stdev_total(), e.stdev_total());
}

TEST(UncertainValueTest, SubtractionFTLM) {
    UncertainValue a(20.0, 0.5, UncertaintySources::FTLM);
    UncertainValue b(15.0, 0.2, UncertaintySources::FTLM);
    UncertainValue c = a - b;
    EXPECT_DOUBLE_EQ(c.mean(), 5.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
    UncertainValue d = a;
    d -= b;
    EXPECT_DOUBLE_EQ(d.mean(), 5.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::FTLM], 0.3);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, SubtractionEXPERIMENT) {
    UncertainValue a(20.0, 0.5, UncertaintySources::EXPERIMENT);
    UncertainValue b(15.0, 0.2, UncertaintySources::EXPERIMENT);
    UncertainValue c = a - b;
    EXPECT_DOUBLE_EQ(c.mean(), 5.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], std::sqrt(0.25+ 0.04));
    UncertainValue d = a;
    d -= b;
    EXPECT_DOUBLE_EQ(d.mean(), 5.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::EXPERIMENT], std::sqrt(0.25+ 0.04));
}

TEST(UncertainValueTest, SubtractionAdditionOfUnaryNegatedFTLM) {
    UncertainValue a(20.0, 0.5, UncertaintySources::FTLM);
    UncertainValue b(15.0, 0.25, UncertaintySources::FTLM);
    UncertainValue c = a - b;
    UncertainValue d = a + (-b);
    EXPECT_DOUBLE_EQ(c.mean(), d.mean());
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], d.stdevs()[UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], d.stdevs()[UncertaintySources::EXPERIMENT]);
    UncertainValue e = a + (-1.0 * b);
    EXPECT_DOUBLE_EQ(c.mean(), e.mean());
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], e.stdevs()[UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], e.stdevs()[UncertaintySources::EXPERIMENT]);
    UncertainValue f = b - a;
    UncertainValue g = b + (-a);
    EXPECT_DOUBLE_EQ(f.mean(), g.mean());
    EXPECT_DOUBLE_EQ(f.stdevs()[UncertaintySources::FTLM], g.stdevs()[UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(f.stdevs()[UncertaintySources::EXPERIMENT], g.stdevs()[UncertaintySources::EXPERIMENT]);
    UncertainValue h = b + (-1.0 * a);
    EXPECT_DOUBLE_EQ(f.mean(), h.mean());
    EXPECT_DOUBLE_EQ(f.stdevs()[UncertaintySources::FTLM], h.stdevs()[UncertaintySources::FTLM]);
    EXPECT_DOUBLE_EQ(f.stdevs()[UncertaintySources::EXPERIMENT], h.stdevs()[UncertaintySources::EXPERIMENT]);
}

TEST(UncertainValueTest, SubtractionCorrelationSignHandlingFTLM) {
    UncertainValue a(20.0, 0.3, UncertaintySources::FTLM);
    UncertainValue b(15.0, 0.5, UncertaintySources::FTLM);
    UncertainValue c(10.0, 0.06, UncertaintySources::FTLM);
    UncertainValue d = (a - b) + c;
    EXPECT_DOUBLE_EQ(d.mean(), 15.0);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::FTLM], 0.14);
    EXPECT_DOUBLE_EQ(d.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
    UncertainValue e = c + (a - b);
    EXPECT_DOUBLE_EQ(e.mean(), 15.0);
    EXPECT_DOUBLE_EQ(e.stdevs()[UncertaintySources::FTLM], 0.14);
    EXPECT_DOUBLE_EQ(e.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
    UncertainValue f = c - (b - a);
    EXPECT_DOUBLE_EQ(f.mean(), 15.0);
    EXPECT_DOUBLE_EQ(f.stdevs()[UncertaintySources::FTLM], 0.14);
    EXPECT_DOUBLE_EQ(f.stdevs()[UncertaintySources::EXPERIMENT], 0.0);

    UncertainValue g = (a - b) + c * 5;
    EXPECT_DOUBLE_EQ(g.mean(), 55.0);
    EXPECT_NEAR(g.stdevs()[UncertaintySources::FTLM], 0.10, 1e-6);
    EXPECT_DOUBLE_EQ(g.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, MultiplicationBasicFTLM) {
    UncertainValue a(3.0, 0.3, UncertaintySources::FTLM);
    UncertainValue b(4.0, 0.4, UncertaintySources::FTLM);
    UncertainValue c = a * b;
    EXPECT_DOUBLE_EQ(c.mean(), 12.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 2.4);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, MultiplicationBasicEXPERIMENT) {
    UncertainValue a(3.0, 0.3, UncertaintySources::EXPERIMENT);
    UncertainValue b(4.0, 0.4, UncertaintySources::EXPERIMENT);
    UncertainValue c = a * b;
    EXPECT_DOUBLE_EQ(c.mean(), 12.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::EXPERIMENT], 12*std::sqrt(0.02));
}

TEST(UncertainValueTest, MultiplicationNegativeFTLM) {
    UncertainValue a(5.0, 0.5, UncertaintySources::FTLM);
    UncertainValue b(-2.0, 0.2, UncertaintySources::FTLM);
    UncertainValue c = a * b;
    EXPECT_DOUBLE_EQ(c.mean(), -10.0);
    EXPECT_DOUBLE_EQ(c.stdevs()[UncertaintySources::FTLM], 2.0);
}

TEST(UncertainValueTest, MultiplicationInverse) {
    UncertainValue a(10.0, {0.3, 0.4});
    UncertainValue b = a / a;
    EXPECT_DOUBLE_EQ(b.mean(), 1.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[UncertaintySources::FTLM], 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[UncertaintySources::EXPERIMENT], std::sqrt(0.0032));
    EXPECT_DOUBLE_EQ(b.stdev_total(), b.stdevs()[UncertaintySources::EXPERIMENT]);
}

TEST(UncertainValueTest, MultiplyByZero) {
    UncertainValue a(10.0, {0.3, 0.4});
    UncertainValue zero(0.0);
    UncertainValue b = a * zero;
    EXPECT_DOUBLE_EQ(b.mean(), 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[FTLM], 0.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[EXPERIMENT], 0.0);
}

TEST(UncertainValueTest, ScalarMultiplication) {
    UncertainValue a(10.0, {0.3, 0.4});
    UncertainValue b = 2.0 * a;
    EXPECT_DOUBLE_EQ(b.mean(), 20.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[FTLM], 0.6);
    EXPECT_DOUBLE_EQ(b.stdevs()[EXPERIMENT], 0.8);
}

TEST(UncertainValueTest, SqrtBasic) {
    UncertainValue a(16.0, 1.6, UncertaintySources::FTLM);
    UncertainValue b = UncertainValue::sqrt(a);
    EXPECT_DOUBLE_EQ(b.mean(), 4.0);
    EXPECT_DOUBLE_EQ(b.stdevs()[UncertaintySources::FTLM], 0.2);
}

TEST(UncertainValueTest, SqrtSquareInversibilityFTLM) {
    UncertainValue a(16.0, 0.25, UncertaintySources::FTLM);
    UncertainValue b = UncertainValue::sqrt(a);
    UncertainValue c = b * b;
    EXPECT_DOUBLE_EQ(a.mean(), c.mean());
    EXPECT_DOUBLE_EQ(a.stdevs()[UncertaintySources::FTLM], c.stdevs()[UncertaintySources::FTLM]);
}

TEST(UncertainValueTest, SquareSqrtInversibilityFTLM) {
    UncertainValue a(16.0, 0.25, UncertaintySources::FTLM);
    UncertainValue b = a * a;
    UncertainValue c = UncertainValue::sqrt(b);
    EXPECT_DOUBLE_EQ(a.mean(), c.mean());
    EXPECT_DOUBLE_EQ(a.stdevs()[UncertaintySources::FTLM], c.stdevs()[UncertaintySources::FTLM]);
}

TEST(UncertainValueTest, InverseIdentity) {
    UncertainValue a(5.0, {0.1, 0.2});
    UncertainValue b = UncertainValue::inv(UncertainValue::inv(a));
    EXPECT_NEAR(a.mean(), b.mean(), 1e-10);
    EXPECT_NEAR(a.stdevs()[FTLM], b.stdevs()[FTLM], 1e-10);
    EXPECT_NEAR(a.stdevs()[EXPERIMENT], b.stdevs()[EXPERIMENT], 1e-10);
}

TEST(UncertainValueTest, MultiplicationSignHandling) {
    UncertainValue pos(5.0, 0.5, UncertaintySources::FTLM);
    UncertainValue neg(-3.0, 0.3, UncertaintySources::FTLM);
    
    UncertainValue pn = pos * neg;
    EXPECT_DOUBLE_EQ(pn.mean(), -15.0);
    EXPECT_NEAR(pn.stdevs()[UncertaintySources::FTLM], 
                3.0, 
                1e-10);

    UncertainValue nn = neg * neg;
    EXPECT_DOUBLE_EQ(nn.mean(), 9.0);
    EXPECT_NEAR(nn.stdevs()[UncertaintySources::FTLM], 
                1.8, 
                1e-10);
}

TEST(UncertainValueTest, DivisionFTLM) {
    UncertainValue a(10.0, 0.5, UncertaintySources::FTLM);
    UncertainValue b(-2.0, 0.1, UncertaintySources::FTLM);

    UncertainValue c = a / b;
    EXPECT_DOUBLE_EQ(c.mean(), -5.0);

    EXPECT_NEAR(c.stdevs()[UncertaintySources::FTLM], 0.0, 1e-10);
}