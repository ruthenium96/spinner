#ifndef SPINNER_ABSTRACTSOLVER_TESTS_H
#define SPINNER_ABSTRACTSOLVER_TESTS_H

#include <random>

#include "gtest/gtest.h"
#include "src/common/runner/Runner.h"
#include "src/nonlinear_solver/AbstractNonlinearSolver.h"

#define RESIDUAL_ERROR_EPSILON 1e-3

template<class T>
std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver> createConcreteSolver();

template<class T>
class fitting_magnetic_susceptibility_simple: public testing::Test {
  protected:
    fitting_magnetic_susceptibility_simple() : solver_(createConcreteSolver<T>()) {}
    std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver> const solver_;
};

TYPED_TEST_SUITE_P(fitting_magnetic_susceptibility_simple);

// TODO: TEST for Theta=0

TYPED_TEST_P(fitting_magnetic_susceptibility_simple, Theta) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> Theta_dist(-100, -1.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double Theta_exact = Theta_dist(rng);
        const double g_exact = 2.0;

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            std::vector<int> mults = {2};
            model::ModelInput model(mults);
            double J_value = Theta_exact;
            auto Theta = model.getSymbols().addSymbol("Theta", J_value);
            model.getSymbols().assignSymbolToTheta(Theta);
            double g_value = g_exact;
            auto g = model.getSymbols().addSymbol("g", g_value, false);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g, i);
            }

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            std::vector<int> mults = {2};
            model::ModelInput model(mults);
            double Theta_value = -10.0;
            auto Theta = model.getSymbols().addSymbol("Theta", Theta_value);
            model.getSymbols().assignSymbolToTheta(Theta);
            double g_value = 2.0;
            auto g = model.getSymbols().addSymbol("g", g_value, false);
            for (size_t i = 0; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g, i);
            }

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(this->solver_);

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double Theta_fitted = runner.getSymbols().getValueOfName(Theta);
            double g_fitted = runner.getSymbols().getValueOfName(g);
            double Theta_range = std::abs(Theta_fitted / 1000);
            double g_range = std::abs(g_fitted / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(Theta_fitted, Theta_exact, Theta_range);
            EXPECT_NEAR(g_fitted, g_exact, g_range);
        }
    }
}

TYPED_TEST_P(fitting_magnetic_susceptibility_simple, fit_theoretical_curve_2222_JAF_g) {
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
            model::ModelInput model(mults);
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

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            std::vector<int> mults = {2, 2, 2, 2};
            model::ModelInput model(mults);
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

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(this->solver_);

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double J_fitted = runner.getSymbols().getValueOfName(J);
            double g_fitted = runner.getSymbols().getValueOfName(g);
            double J_range = std::abs(J_fitted / 1000);
            double g_range = std::abs(g_fitted / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(J_fitted, J_exact, J_range);
            EXPECT_NEAR(g_fitted, g_exact, g_range);
        }
    }
}

TYPED_TEST_P(fitting_magnetic_susceptibility_simple, fit_theoretical_curve_222222_JAF_g) {
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
            model::ModelInput model(mults);
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

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            std::vector<int> mults = {2, 2, 2, 2, 2, 2};
            model::ModelInput model(mults);
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

            runner::Runner runner(model);

            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(this->solver_);

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double J_fitted = runner.getSymbols().getValueOfName(J);
            double g_fitted = runner.getSymbols().getValueOfName(g);
            double J_range = std::abs(J_fitted / 1000);
            double g_range = std::abs(g_fitted / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(J_fitted, J_exact, J_range);
            EXPECT_NEAR(g_fitted, g_exact, g_range);
        }
    }
}

TYPED_TEST_P(fitting_magnetic_susceptibility_simple, fit_theoretical_curve_222222_JFM_g) {
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
            model::ModelInput model(mults);
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

            runner::Runner runner(model);
            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            std::vector<int> mults = {2, 2, 2, 2, 2, 2};
            model::ModelInput model(mults);
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

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(this->solver_);

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double J_fitted = runner.getSymbols().getValueOfName(J);
            double g_fitted = runner.getSymbols().getValueOfName(g);
            double J_range = std::abs(J_fitted / 1000);
            double g_range = std::abs(g_fitted / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(J_fitted, J_exact, J_range);
            EXPECT_NEAR(g_fitted, g_exact, g_range);
        }
    }
}

TYPED_TEST_P(fitting_magnetic_susceptibility_simple, fit_theoretical_curve_222222_g_fixed_g) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> g_one_dist(1.8, 2.0);
    std::uniform_real_distribution<double> g_two_dist(2.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double g_one_exact = g_one_dist(rng);
        const double g_two_exact = g_two_dist(rng);
        std::vector<int> mults = {2, 2, 2, 2, 2, 2};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            model::ModelInput model(mults);
            double g_one_value = g_one_exact;
            double g_two_value = g_two_exact;
            auto g_one = model.getSymbols().addSymbol("g1", g_one_value);
            auto g_two = model.getSymbols().addSymbol("g2", g_two_value);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.getSymbols().assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            model::ModelInput model(mults);
            double g_two_value = 3.0;
            auto g_one = model.getSymbols().addSymbol("g1", g_one_exact, false);
            auto g_two = model.getSymbols().addSymbol("g2", g_two_value);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.getSymbols().assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(this->solver_);

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double g_one_fitted = runner.getSymbols().getValueOfName(g_one);
            double g_two_fitted = runner.getSymbols().getValueOfName(g_two);
            double g_one_range = std::abs((g_one_fitted) / 1000);
            double g_two_range = std::abs((g_two_fitted) / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(g_one_fitted, g_one_exact, g_one_range);
            EXPECT_NEAR(g_two_fitted, g_two_exact, g_two_range);
        }
    }
}

TYPED_TEST_P(fitting_magnetic_susceptibility_simple, fit_theoretical_curve_222222_JAF_g_fixed_g) {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> J_dist(-200, -0.1);
    std::uniform_real_distribution<double> g_one_dist(1.8, 2.0);
    std::uniform_real_distribution<double> g_two_dist(2.8, 3.0);

    for (size_t _ = 0; _ < 20; ++_) {
        const double J_exact = J_dist(rng);
        const double g_one_exact = g_one_dist(rng);
        const double g_two_exact = g_two_dist(rng);
        std::vector<int> mults = {2, 2, 2, 2, 2, 2};

        std::vector<magnetic_susceptibility::ValueAtTemperature> values;

        {
            model::ModelInput model(mults);
            double J_value = J_exact;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 4)
                .assignSymbolToIsotropicExchange(J, 4, 5)
                .assignSymbolToIsotropicExchange(J, 5, 0);
            double g_one_value = g_one_exact;
            double g_two_value = g_two_exact;
            auto g_one = model.getSymbols().addSymbol("g1", g_one_value);
            auto g_two = model.getSymbols().addSymbol("g2", g_two_value);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.getSymbols().assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);

            runner.BuildSpectra();
            runner.BuildMuSquaredWorker();

            for (size_t i = 1; i < 301; ++i) {
                magnetic_susceptibility::ValueAtTemperature value_at_temperature = {
                    static_cast<double>(i),
                    runner.getMagneticSusceptibilityController().calculateTheoreticalMuSquared(i)};
                values.push_back(value_at_temperature);
            }
        }

        {
            model::ModelInput model(mults);
            double J_value = -10.0;
            auto J = model.getSymbols().addSymbol("J", J_value);
            model.getSymbols()
                .assignSymbolToIsotropicExchange(J, 0, 1)
                .assignSymbolToIsotropicExchange(J, 1, 2)
                .assignSymbolToIsotropicExchange(J, 2, 3)
                .assignSymbolToIsotropicExchange(J, 3, 4)
                .assignSymbolToIsotropicExchange(J, 4, 5)
                .assignSymbolToIsotropicExchange(J, 5, 0);
            double g_two_value = 3.0;
            auto g_one = model.getSymbols().addSymbol("g1", g_one_exact, false);
            auto g_two = model.getSymbols().addSymbol("g2", g_two_value);
            for (size_t i = 0; i < mults.size() / 2; ++i) {
                model.getSymbols().assignSymbolToGFactor(g_one, i);
            }
            for (size_t i = mults.size() / 2; i < mults.size(); ++i) {
                model.getSymbols().assignSymbolToGFactor(g_two, i);
            }

            runner::Runner runner(model);
            runner.initializeExperimentalValues(
                values,
                magnetic_susceptibility::mu_squared_in_bohr_magnetons_squared,
                1);

            runner.minimizeResidualError(this->solver_);

            double residual_error =
                runner.getMagneticSusceptibilityController().calculateResidualError();
            double g_one_fitted = runner.getSymbols().getValueOfName(g_one);
            double g_two_fitted = runner.getSymbols().getValueOfName(g_two);
            double g_one_range = std::abs((g_one_fitted) / 1000);
            double g_two_range = std::abs((g_two_fitted) / 1000);

            EXPECT_NEAR(residual_error, 0, RESIDUAL_ERROR_EPSILON);
            EXPECT_NEAR(g_one_fitted, g_one_exact, g_one_range);
            EXPECT_NEAR(g_two_fitted, g_two_exact, g_two_range);
        }
    }
}

REGISTER_TYPED_TEST_SUITE_P(
    fitting_magnetic_susceptibility_simple,
    Theta,
    fit_theoretical_curve_2222_JAF_g,
    fit_theoretical_curve_222222_JFM_g,
    fit_theoretical_curve_222222_JAF_g,
    fit_theoretical_curve_222222_g_fixed_g,
    fit_theoretical_curve_222222_JAF_g_fixed_g);

#endif  //SPINNER_ABSTRACTSOLVER_TESTS_H
