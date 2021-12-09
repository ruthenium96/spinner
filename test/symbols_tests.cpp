#include "common/runner/Runner.h"
#include "gtest/gtest.h"

TEST(symbols, throw_2222_isotropic_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J_value = 10;
    auto tripledJ = runner.AddSymbol("3J", 3 * J_value);
    auto J = runner.AddSymbol("J", J_value);
    auto doubledJ = runner.AddSymbol("2J", 2 * J_value);
    runner.AddIsotropicExchange(tripledJ, 0, 1);
    runner.AddIsotropicExchange(J, 1, 2);
    runner.AddIsotropicExchange(doubledJ, 2, 3);
    runner.AddIsotropicExchange(J, 3, 0);
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
    auto firstJ = runner.AddSymbol("J1", J);
    auto secondJ = runner.AddSymbol("J2", J);
    runner.AddIsotropicExchange(firstJ, 0, 1);
    runner.AddIsotropicExchange(firstJ, 1, 2);
    runner.AddIsotropicExchange(firstJ, 2, 3);
    runner.AddIsotropicExchange(secondJ, 3, 0);
    runner.FinalizeIsotropicInteraction();

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(symbols, throw_2222_gfactor_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J_value = 10;
    auto J = runner.AddSymbol("J", J_value);
    auto firstg = runner.AddSymbol("g1", 2.0);
    auto secondg = runner.AddSymbol("g2", 3.0);
    runner.AddIsotropicExchange(J, 0, 1);
    runner.AddIsotropicExchange(J, 1, 2);
    runner.AddIsotropicExchange(J, 2, 3);
    runner.AddIsotropicExchange(J, 3, 0);
    runner.AddGFactor(firstg, 0);
    runner.AddGFactor(firstg, 1);
    runner.AddGFactor(firstg, 2);
    runner.AddGFactor(secondg, 3);

    runner.FinalizeIsotropicInteraction();

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(symbols, throw_2222_gfactor_accidental_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J_value = 10;
    double g_value = 2.0;
    auto J = runner.AddSymbol("J", J_value);
    auto firstg = runner.AddSymbol("g1", g_value);
    auto secondg = runner.AddSymbol("g2", g_value);
    runner.AddIsotropicExchange(J, 0, 1);
    runner.AddIsotropicExchange(J, 1, 2);
    runner.AddIsotropicExchange(J, 2, 3);
    runner.AddIsotropicExchange(J, 3, 0);
    runner.AddGFactor(firstg, 0);
    runner.AddGFactor(firstg, 1);
    runner.AddGFactor(firstg, 2);
    runner.AddGFactor(secondg, 3);

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

TEST(symbols, throw_set_new_value_to_unchangeable_symbol) {
    symbols::Symbols symbols_(10);

    auto unchangeable = symbols_.addSymbol("Unchangeable", NAN, false, symbols::not_specified);

    EXPECT_THROW(
        symbols_.setNewValueToChangeableSymbol(unchangeable, INFINITY),
        std::invalid_argument);
}

TEST(symbols, throw_specified_not_as_J_symbol) {
    size_t number_of_spins = 10;
    symbols::Symbols symbols_(number_of_spins);

    double gChangeable = 2;

    auto not_specified = symbols_.addSymbol("not_specified", NAN, true, symbols::not_specified);
    symbols_.addIsotropicExchange(not_specified, 1, 3);
    auto g_changeable = symbols_.addSymbol("gChangeable", gChangeable, true, symbols::g_factor);
    EXPECT_THROW(symbols_.addIsotropicExchange(g_changeable, 2, 7), std::invalid_argument);
}

TEST(symbols, throw_specified_not_as_g_symbol) {
    size_t number_of_spins = 10;
    symbols::Symbols symbols_(number_of_spins);

    double JChangeable = 10;

    auto not_specified = symbols_.addSymbol("not_specified", NAN, true, symbols::not_specified);
    symbols_.addGFactor(not_specified, 1);
    auto J_changeable = symbols_.addSymbol("JChangeable", JChangeable, true, symbols::J);
    EXPECT_THROW(symbols_.addGFactor(J_changeable, 2), std::invalid_argument);
}

TEST(symbols, set_new_value_to_changeable_J_g) {
    size_t number_of_spins = 10;
    symbols::Symbols symbols_(number_of_spins);

    double JChangeable_value = 10;
    double gChangeable_value = 2;

    auto J_changeable = symbols_.addSymbol("JChangeable", JChangeable_value, true, symbols::J);
    auto g_changeable =
        symbols_.addSymbol("gChangeable", gChangeable_value, true, symbols::g_factor);
    symbols_.addIsotropicExchange(J_changeable, 2, 7);
    for (size_t i = 0; i < number_of_spins; ++i) {
        symbols_.addGFactor(g_changeable, i);
    }
    auto shared_ptr_J = symbols_.constructIsotropicExchangeParameters();
    auto shared_ptr_g = symbols_.constructGFactorParameters();
    EXPECT_EQ(JChangeable_value, shared_ptr_J->operator()(2, 7));
    EXPECT_EQ(gChangeable_value, shared_ptr_g->operator()(7));
    symbols_.setNewValueToChangeableSymbol(J_changeable, 2 * JChangeable_value);
    EXPECT_EQ(JChangeable_value, shared_ptr_J->operator()(2, 7));
    symbols_.setNewValueToChangeableSymbol(g_changeable, 2 * gChangeable_value);
    EXPECT_EQ(gChangeable_value, shared_ptr_g->operator()(7));
    symbols_.updateIsotropicExchangeParameters();
    EXPECT_EQ(2 * JChangeable_value, shared_ptr_J->operator()(2, 7));
    symbols_.updateGFactorParameters();
    EXPECT_EQ(2 * gChangeable_value, shared_ptr_g->operator()(7));
}