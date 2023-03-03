#ifndef SPINNER_3D_NIT2_H
#define SPINNER_3D_NIT2_H

#include "src/common/runner/Runner.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"

// Calculation of chiT of a simple diradical from https://doi.org/10.1039/B414284E
void three_d_nit_two() {
    // Values: cm^3K/mol. Data was extracted from Figure 4 of the article.
    std::vector<magnetic_susceptibility::ValueAtTemperature> values = {
        {295.135, 0.616}, {284.865, 0.614}, {275.135, 0.608}, {47.568, 0.186},  {55.676, 0.229},
        {65.405, 0.284},  {75.135, 0.330},  {84.865, 0.372},  {94.595, 0.406},  {104.324, 0.433},
        {115.135, 0.456}, {125.405, 0.480}, {135.135, 0.499}, {43.243, 0.158},  {40.541, 0.129},
        {36.216, 0.103},  {31.892, 0.074},  {27.568, 0.050},  {23.784, 0.029},  {265.405, 0.607},
        {255.676, 0.603}, {244.865, 0.598}, {234.595, 0.594}, {224.865, 0.588}, {215.135, 0.583},
        {205.405, 0.579}, {195.135, 0.569}, {185.405, 0.561}, {174.595, 0.552}, {164.865, 0.541},
        {154.595, 0.527}, {144.865, 0.515}};
    // There are only two paramagnetic centers.
    std::vector<int> mults = {2, 2};
    model::ModelInput model(mults);
    // Start from reasonable guess:
    double J_initial = -50.0;
    auto J = model.modifySymbolicWorker().addSymbol("J", J_initial);
    model.modifySymbolicWorker().assignSymbolToIsotropicExchange(J, 0, 1);
    double g_initial = 2.0;
    auto g = model.modifySymbolicWorker().addSymbol("g", g_initial);
    model.modifySymbolicWorker().assignSymbolToGFactor(g, 0);
    model.modifySymbolicWorker().assignSymbolToGFactor(g, 1);
    double theta_initial = -10;
    auto theta = model.modifySymbolicWorker().addSymbol("theta", theta_initial);
    model.modifySymbolicWorker().assignSymbolToTheta(theta);

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

#endif  //SPINNER_3D_NIT2_H
