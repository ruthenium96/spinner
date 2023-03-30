#ifndef SPINNER_DTD_NIT2_H
#define SPINNER_DTD_NIT2_H

#include "src/common/runner/Runner.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"

// Calculation of chiT of a simple diradical from https://doi.org/10.1039/B414284E
void dtd_nit_two() {
    // Values: cm^3K/mol. Data was extracted from Figure 7 of the article.
    std::vector<magnetic_susceptibility::ValueAtTemperature> values = {
        {289.868, 0.751}, {280.161, 0.750}, {269.973, 0.746}, {259.295, 0.744}, {250.079, 0.738},
        {239.890, 0.734}, {229.214, 0.731}, {219.994, 0.728}, {209.321, 0.722}, {199.133, 0.718},
        {179.240, 0.709}, {159.348, 0.700}, {138.987, 0.678}, {119.116, 0.651}, {99.252, 0.619},
        {78.446, 0.563},  {68.782, 0.526},  {58.643, 0.481},  {48.044, 0.415},  {38.930, 0.325},
        {28.880, 0.206},  {27.454, 0.181},  {24.577, 0.151},  {22.667, 0.125},  {20.758, 0.099},
        {18.361, 0.074},  {15.004, 0.040},
    };
    // There are only two paramagnetic centers.
    std::vector<spin_algebra::Multiplicity> mults = {2, 2};
    model::ModelInput model(mults);
    // Start from reasonable guess:
    double J_initial = -50.0;
    auto J = model.modifySymbolicWorker().addSymbol("J", J_initial);
    model.modifySymbolicWorker().assignSymbolToIsotropicExchange(J, 0, 1);
    double g_initial = 2.0;
    auto g = model.modifySymbolicWorker().addSymbol("g", g_initial);
    model.modifySymbolicWorker().assignSymbolToGFactor(g, 0);
    model.modifySymbolicWorker().assignSymbolToGFactor(g, 1);

    runner::Runner runner(model);
    runner.initializeExperimentalValues(
        values,
        magnetic_susceptibility::chiT_in_cm_cubed_kelvin_per_mol,
        1);

    runner.minimizeResidualError(std::make_shared<nonlinear_solver::optimNMAdapter>());

    auto theor_final = runner.getMagneticSusceptibilityController().getTheoreticalValues();

    std::cout << std::endl;
    for (auto [temperature, value] : theor_final) {
        std::cout << temperature << " " << value << std::endl;
    }

    // The answer of calculation:
    //
    // J: -48.2693
    // g: 1.98555
    // theta: -24.7272
    // RSS = 0.034688

    // Fit of chiT from the article: J = -47,5 K, Theta = -25 K.
    // ESR data: g = 2.0065, J = -43 K
}

#endif  //SPINNER_DTD_NIT2_H
