#include <random>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

TEST(magnetic_susceptibility, fit_theoretical_curve_222222_JAF_g) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-200, -0.1);
    std::uniform_real_distribution<double> g_dist(1.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double g_exact = g_dist(rng);

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            std::vector<int> mults = {2, 2, 2, 2, 2, 2};
            model::Model model(mults);
            double J_value = J_exact;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 4)
                .assignSymbolToIsotropicExchange(J, 4, 5)
                .assignSymbolToIsotropicExchange(J, 5, 0);

            double g_value = g_exact;
            auto g = model.getSymbols().addSymbol("g", g_value);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g, i);
            }

            model.InitializeSSquared();
            model.InitializeIsotropicExchangeDerivatives();

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMuSquaredWorker().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            std::vector<int> mults = {2, 2, 2, 2, 2, 2};
            model::Model model(mults);
            double J_value = -10.0;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 4)
                .assignSymbolToIsotropicExchange(J, 4, 5)
                .assignSymbolToIsotropicExchange(J, 5, 0);

            double g_value = g_exact;
            auto g = model.getSymbols().addSymbol("g", g_value);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g, i);
            }

            model.InitializeSSquared();
            model.InitializeIsotropicExchangeDerivatives();

            runner::Runner runner(model);

            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError();

            double residual_error = runner.getMuSquaredWorker().calculateResidualError();
            double J_fitted = runner.getSymbols().getValueOfName(J);
            double g_fitted = runner.getSymbols().getValueOfName(g);
            double J_range = std::abs(J_fitted / 1000);
            double g_range = std::abs(g_fitted / 1000);

            // TODO: epsilon
            EXPECT_NEAR(residual_error, 0, 1e-3);
            EXPECT_NEAR(J_fitted, J_exact, J_range);
            EXPECT_NEAR(g_fitted, g_exact, g_range);
        }
    }
}

TEST(magnetic_susceptibility, fit_theoretical_curve_222222_JFM_g) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(20, 200);
    std::uniform_real_distribution<double> g_dist(1.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double g_exact = g_dist(rng);

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            std::vector<int> mults = {2, 2, 2, 2, 2, 2};
            model::Model model(mults);
            double J_value = J_exact;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 4)
                .assignSymbolToIsotropicExchange(J, 4, 5)
                .assignSymbolToIsotropicExchange(J, 5, 0);

            double g_value = g_exact;
            auto g = model.getSymbols().addSymbol("g", g_value);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g, i);
            }
            model.InitializeSSquared();
            model.InitializeIsotropicExchangeDerivatives();

            runner::Runner runner(model);
            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMuSquaredWorker().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            std::vector<int> mults = {2, 2, 2, 2, 2, 2};
            model::Model model(mults);
            double J_value = 40;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 4)
                .assignSymbolToIsotropicExchange(J, 4, 5)
                .assignSymbolToIsotropicExchange(J, 5, 0);

            double g_value = g_exact;
            auto g = model.getSymbols().addSymbol("g", g_value);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g, i);
            }
            model.InitializeSSquared();
            model.InitializeIsotropicExchangeDerivatives();

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError();

            double residual_error = runner.getMuSquaredWorker().calculateResidualError();
            double J_fitted = runner.getSymbols().getValueOfName(J);
            double g_fitted = runner.getSymbols().getValueOfName(g);
            double J_range = std::abs(J_fitted / 1000);
            double g_range = std::abs(g_fitted / 1000);

            // TODO: epsilon
            EXPECT_NEAR(residual_error, 0, 1e-3);
            EXPECT_NEAR(J_fitted, J_exact, J_range);
            EXPECT_NEAR(g_fitted, g_exact, g_range);
        }
    }
}