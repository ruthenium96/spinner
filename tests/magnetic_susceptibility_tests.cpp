#include <random>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"

TEST(magnetic_susceptibility, throw_experimental_values_worker_set_values_length_error) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    std::vector<magnetic_susceptibility::ValueAtTemperature> theor_values = exp_values;
    theor_values.pop_back();

    EXPECT_THROW(
        experimental_values_worker.setTheoreticalMuSquared(theor_values),
        std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_get_values_empty) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    EXPECT_THROW(experimental_values_worker.getTheoreticalValues(), std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_residual_error_empty) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    EXPECT_THROW(experimental_values_worker.calculateResidualError(), std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_derivative_empty) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 1}, {2, 2}, {3, 3}};

    magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    EXPECT_THROW(experimental_values_worker.calculateDerivative(), std::length_error);
}

TEST(magnetic_susceptibility, throw_experimental_values_worker_empty_experimental_values) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {};

    EXPECT_THROW(
        magnetic_susceptibility::ExperimentalValuesWorker experimental_values_worker(
            exp_values,
            magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
            1),
        std::length_error);
}

TEST(magnetic_susceptibility, do_not_throw_experimental_before_theoretical) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 20}};
    std::vector<int> mults = {2, 2};
    model::Model model(mults);
    double J_value = 10;
    auto J = model.getSymbols().addSymbol("J", J_value);
    model.getSymbols().assignSymbolToIsotropicExchange(J, 0, 1);
    double g_value = 2.0;
    auto g = model.getSymbols().addSymbol("g", g_value);
    for (size_t i = 0; i < mults.size(); ++i) {
        model.getSymbols().assignSymbolToGFactor(g, i);
    }
    model.InitializeSSquared();
    model.InitializeIsotropicExchangeDerivatives();

    runner::Runner runner(model);

    runner.initializeExperimentalValues(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1);

    runner.BuildSpectra();
    EXPECT_NO_THROW(runner.BuildMuSquaredWorker());
}

TEST(magnetic_susceptibility, do_not_throw_theoretical_before_experimental) {
    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values = {{1, 20}};
    std::vector<int> mults = {2, 2};
    model::Model model(mults);
    double J_value = 10;
    auto J = model.getSymbols().addSymbol("J", J_value);
    model.getSymbols().assignSymbolToIsotropicExchange(J, 0, 1);
    double g_value = 2.0;
    auto g = model.getSymbols().addSymbol("g", g_value);
    for (size_t i = 0; i < mults.size(); ++i) {
        model.getSymbols().assignSymbolToGFactor(g, i);
    }
    model.InitializeSSquared();
    model.InitializeIsotropicExchangeDerivatives();

    runner::Runner runner(model);

    runner.BuildSpectra();
    runner.BuildMuSquaredWorker();
    EXPECT_NO_THROW(runner.initializeExperimentalValues(
        exp_values,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        1));
}

TEST(magnetic_susceptibility, value_mu_squared_reversibility) {
    std::vector<magnetic_susceptibility::ExperimentalValuesEnum> values_enum = {
        magnetic_susceptibility::mu_in_bohr_magnetons,
        magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
        magnetic_susceptibility::chiT_in_cm_cubed_kelvin_per_mol};

    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist(1, 100);

    std::vector<magnetic_susceptibility::ValueAtTemperature> exp_values;

    for (size_t i = 0; i < 300; ++i) {
        magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
            static_cast<double>(i),
            static_cast<double>(dist(rng))};
        exp_values.push_back(value_at_temperature);
    }

    for (const auto& value_enum : values_enum) {
        for (size_t ratio = 1; ratio < 10; ++ratio) {
            magnetic_susceptibility::ExperimentalValuesWorker worker(exp_values, value_enum, ratio);
            worker.setTheoreticalMuSquared(worker.getExperimentalMuSquared());
            auto double_transformed = worker.getTheoreticalValues();
            for (size_t i = 0; i < double_transformed.size(); ++i) {
                EXPECT_DOUBLE_EQ(double_transformed[i].value, exp_values[i].value);
            }
        }
    }
}

TEST(magnetic_susceptibility, fit_theoretical_curve_2222_JAF_g) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-200, -0.1);
    std::uniform_real_distribution<double> g_dist(1.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double g_exact = g_dist(rng);

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            std::vector<int> mults = {2, 2, 2, 2};
            model::Model model(mults);
            double J_value = J_exact;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);
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
            std::vector<int> mults = {2, 2, 2, 2};
            model::Model model(mults);
            double J_value = -10.0;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 0);
            double g_value = 2.0;
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