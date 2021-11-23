#include "common/runner/Runner.h"
#include "gtest/gtest.h"

TEST(symbols, throw_2222_isotropic_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J = 10;
    runner.AddSymbol("3J", 3 * J);
    runner.AddSymbol("J", J);
    runner.AddSymbol("2J", 2 * J);
    runner.AddIsotropicExchange("3J", 0, 1);
    runner.AddIsotropicExchange("J", 1, 2);
    runner.AddIsotropicExchange("2J", 2, 3);
    runner.AddIsotropicExchange("J", 3, 0);
    runner.FinalizeIsotropicInteraction();

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(symbols, throw_2222_isotropic_accidental_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J = 10;
    runner.AddSymbol("J1", J);
    runner.AddSymbol("J2", J);
    runner.AddIsotropicExchange("J1", 0, 1);
    runner.AddIsotropicExchange("J1", 1, 2);
    runner.AddIsotropicExchange("J1", 2, 3);
    runner.AddIsotropicExchange("J2", 3, 0);
    runner.FinalizeIsotropicInteraction();

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(symbols, throw_2222_gfactor_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J = 10;
    runner.AddSymbol("J", J);
    runner.AddSymbol("g1", 2.0);
    runner.AddSymbol("g2", 3.0);
    runner.AddIsotropicExchange("J", 0, 1);
    runner.AddIsotropicExchange("J", 1, 2);
    runner.AddIsotropicExchange("J", 2, 3);
    runner.AddIsotropicExchange("J", 3, 0);
    runner.AddGFactor("g1", 0);
    runner.AddGFactor("g1", 1);
    runner.AddGFactor("g1", 2);
    runner.AddGFactor("g2", 3);

    runner.FinalizeIsotropicInteraction();

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(symbols, throw_2222_gfactor_accidental_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J = 10;
    double g = 2.0;
    runner.AddSymbol("J", J);
    runner.AddSymbol("g1", g);
    runner.AddSymbol("g2", g);
    runner.AddIsotropicExchange("J", 0, 1);
    runner.AddIsotropicExchange("J", 1, 2);
    runner.AddIsotropicExchange("J", 2, 3);
    runner.AddIsotropicExchange("J", 3, 0);
    runner.AddGFactor("g1", 0);
    runner.AddGFactor("g1", 1);
    runner.AddGFactor("g1", 2);
    runner.AddGFactor("g2", 3);

    runner.FinalizeIsotropicInteraction();

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

// TODO: make these tests pass
//TEST(symbols, throw_2222_gfactor_all_were_not_initialized) {
//    std::vector<int> mults = {2, 2, 2, 2};
//
//    runner::Runner runner(mults);
//
//    runner.TzSort();
//
//    double J = 10;
//    double g = 2.0;
//    runner.AddSymbol("J", J);
//    runner.AddSymbol("g1", g);
//    runner.AddIsotropicExchange("J", 0, 1);
//    runner.AddIsotropicExchange("J", 1, 2);
//    runner.AddIsotropicExchange("J", 2, 3);
//    runner.AddIsotropicExchange("J", 3, 0);
//
//    runner.FinalizeIsotropicInteraction();
//
//    EXPECT_THROW(runner.BuildMatrices(), std::length_error);
//}
//
//TEST(symbols, throw_2222_gfactor_any_was_not_initialized) {
//    std::vector<int> mults = {2, 2, 2, 2};
//
//    runner::Runner runner(mults);
//
//    runner.TzSort();
//
//    double J = 10;
//    double g = 2.0;
//    runner.AddSymbol("J", J);
//    runner.AddSymbol("g1", g);
//    runner.AddIsotropicExchange("J", 0, 1);
//    runner.AddIsotropicExchange("J", 1, 2);
//    runner.AddIsotropicExchange("J", 2, 3);
//    runner.AddIsotropicExchange("J", 3, 0);
//    runner.AddGFactor("g1", 0);
//    runner.AddGFactor("g1", 1);
//    runner.AddGFactor("g1", 2);
//
//    runner.FinalizeIsotropicInteraction();
//
//    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
//}
