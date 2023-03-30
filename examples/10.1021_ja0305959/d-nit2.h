#ifndef SPINNER_D_NIT2_H
#define SPINNER_D_NIT2_H

#include "src/common/runner/Runner.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"

// Calculation of chiT of a simple diradical from https://doi.org/10.1021/ja0305959
void d_nit_two() {
    // Values: cm^3K/mol. Data was extracted from Figure 2a of the article.
    std::vector<magnetic_susceptibility::ValueAtTemperature> values = {
        {288.445, 0.372}, {278.239, 0.360}, {269.433, 0.347}, {258.300, 0.334}, {248.569, 0.318},
        {238.378, 0.300}, {229.583, 0.283}, {219.395, 0.263}, {209.206, 0.244}, {199.484, 0.225},
        {189.760, 0.206}, {179.111, 0.185}, {169.856, 0.165}, {159.678, 0.142}, {149.961, 0.120},
        {140.239, 0.100}, {130.053, 0.080}, {120.331, 0.061}, {109.671, 0.044}, {99.935, 0.031},
        {89.731, 0.019},  {79.981, 0.011},  {70.692, 0.006},  {60.464, 0.003},  {46.046, 0.004},
        {35.349, 0.003},  {24.186, 0.003},  {13.023, 0.003},  {5.581, 0.003}};
    // There are only two paramagnetic centers.
    std::vector<spin_algebra::Multiplicity> mults = {2, 2};
    model::ModelInput model(mults);
    // Start from reasonable guess:
    double J_initial = -200.0;
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
    // J: -229.865
    // g: 1.98444
    // RSS = 0.00658317

    // It should be mentioned, that Spinner and the article use different form of
    // isotropic Hamiltonian, so J's from the article should be divided by two.
    // Fit of chiT from the article: J = -230 K
    // ESR data: g = 2.007, J = -240 K
}

#endif  //SPINNER_D_NIT2_H
