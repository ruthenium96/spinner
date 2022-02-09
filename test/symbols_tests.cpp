#include "common/runner/Runner.h"
#include "gtest/gtest.h"

TEST(symbols, throw_add_the_same_symbol_name) {
    symbols::Symbols symbols(2);
    symbols.addSymbol("same", 10);
    EXPECT_THROW(symbols.addSymbol("same", 10), std::invalid_argument);
}

TEST(symbols, throw_assign_the_same_exchange) {
    symbols::Symbols symbols(2);
    auto J = symbols.addSymbol("same", 10, symbols::J);
    symbols.assignSymbolToIsotropicExchange(J, 0, 1);
    EXPECT_THROW(symbols.assignSymbolToIsotropicExchange(J, 0, 1), std::invalid_argument);
}

TEST(symbols, throw_2222_isotropic_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J_value = 10;
    auto tripledJ = runner.getMutableSymbols().addSymbol("3J", 3 * J_value);
    auto J = runner.getMutableSymbols().addSymbol("J", J_value);
    auto doubledJ = runner.getMutableSymbols().addSymbol("2J", 2 * J_value);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(tripledJ, 0, 1);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 1, 2);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(doubledJ, 2, 3);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 3, 0);

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(symbols, throw_2222_isotropic_accidental_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J = 10;
    auto firstJ = runner.getMutableSymbols().addSymbol("J1", J);
    auto secondJ = runner.getMutableSymbols().addSymbol("J2", J);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(firstJ, 0, 1);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(firstJ, 1, 2);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(firstJ, 2, 3);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(secondJ, 3, 0);

    EXPECT_THROW(runner.BuildMatrices(), std::invalid_argument);
}

TEST(symbols, throw_2222_gfactor_inconsistent_symmetry) {
    std::vector<int> mults = {2, 2, 2, 2};

    runner::Runner runner(mults);

    runner.TzSort();
    runner.Symmetrize(group::Group::S2, {{1, 0, 3, 2}});
    runner.Symmetrize(group::Group::S2, {{3, 2, 1, 0}});

    double J_value = 10;
    auto J = runner.getMutableSymbols().addSymbol("J", J_value);
    auto firstg = runner.getMutableSymbols().addSymbol("g1", 2.0);
    auto secondg = runner.getMutableSymbols().addSymbol("g2", 3.0);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 0, 1);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 1, 2);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 2, 3);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 3, 0);
    runner.getMutableSymbols().assignSymbolToGFactor(firstg, 0);
    runner.getMutableSymbols().assignSymbolToGFactor(firstg, 1);
    runner.getMutableSymbols().assignSymbolToGFactor(firstg, 2);
    runner.getMutableSymbols().assignSymbolToGFactor(secondg, 3);

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
    auto J = runner.getMutableSymbols().addSymbol("J", J_value);
    auto firstg = runner.getMutableSymbols().addSymbol("g1", g_value);
    auto secondg = runner.getMutableSymbols().addSymbol("g2", g_value);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 0, 1);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 1, 2);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 2, 3);
    runner.getMutableSymbols().assignSymbolToIsotropicExchange(J, 3, 0);
    runner.getMutableSymbols().assignSymbolToGFactor(firstg, 0);
    runner.getMutableSymbols().assignSymbolToGFactor(firstg, 1);
    runner.getMutableSymbols().assignSymbolToGFactor(firstg, 2);
    runner.getMutableSymbols().assignSymbolToGFactor(secondg, 3);

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
//    runner.modifySymbol().addSymbol("J", J);
//    runner.modifySymbol().addSymbol("g1", g);
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
//    runner.modifySymbol().addSymbol("J", J);
//    runner.modifySymbol().addSymbol("g1", g);
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
    symbols_.assignSymbolToIsotropicExchange(not_specified, 1, 3);
    auto g_changeable = symbols_.addSymbol("gChangeable", gChangeable, true, symbols::g_factor);
    EXPECT_THROW(
        symbols_.assignSymbolToIsotropicExchange(g_changeable, 2, 7),
        std::invalid_argument);
}

TEST(symbols, throw_specified_not_as_g_symbol) {
    size_t number_of_spins = 10;
    symbols::Symbols symbols_(number_of_spins);

    double JChangeable = 10;

    auto not_specified = symbols_.addSymbol("not_specified", NAN, true, symbols::not_specified);
    symbols_.assignSymbolToGFactor(not_specified, 1);
    auto J_changeable = symbols_.addSymbol("JChangeable", JChangeable, true, symbols::J);
    EXPECT_THROW(symbols_.assignSymbolToGFactor(J_changeable, 2), std::invalid_argument);
}

TEST(symbols, set_new_value_to_changeable_J_g) {
    size_t number_of_spins = 10;
    symbols::Symbols symbols_(number_of_spins);

    double JChangeable_value = 10;
    double gChangeable_value = 2;

    auto J_changeable = symbols_.addSymbol("JChangeable", JChangeable_value, true, symbols::J);
    auto g_changeable =
        symbols_.addSymbol("gChangeable", gChangeable_value, true, symbols::g_factor);
    symbols_.assignSymbolToIsotropicExchange(J_changeable, 2, 7);
    for (size_t i = 0; i < number_of_spins; ++i) {
        symbols_.assignSymbolToGFactor(g_changeable, i);
    }
    auto shared_ptr_J = symbols_.getIsotropicExchangeParameters();
    auto shared_ptr_g = symbols_.getGFactorParameters();
    EXPECT_EQ(JChangeable_value, shared_ptr_J->operator()(2, 7));
    EXPECT_EQ(gChangeable_value, shared_ptr_g->operator()(7));
    symbols_.setNewValueToChangeableSymbol(J_changeable, 2 * JChangeable_value);
    EXPECT_EQ(2 * JChangeable_value, shared_ptr_J->operator()(2, 7));
    EXPECT_EQ(gChangeable_value, shared_ptr_g->operator()(7));
    symbols_.setNewValueToChangeableSymbol(g_changeable, 2 * gChangeable_value);
    EXPECT_EQ(2 * JChangeable_value, shared_ptr_J->operator()(2, 7));
    EXPECT_EQ(2 * gChangeable_value, shared_ptr_g->operator()(7));
}