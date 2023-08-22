#include <utility>

#include "src/common/runner/Runner.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"
#include "src/nonlinear_solver/stlbfgs/stlbfgsAdapter.h"

// Calculation of chiT of nickel crystal from https://doi.org/10.1039/D0DT03184D

std::vector<magnetic_susceptibility::ValueAtTemperature> experimental_values =  // cm^3K/mol
    {{1.99979, 1.055030643},    {1.99946, 1.054682046625}, {2.49832, 1.1990345225},
     {2.99917, 1.30990060875},  {3.49925, 1.3940575775},   {4.01245, 1.46086179125},
     {4.49819, 1.52034946875},  {4.9995, 1.571011065},     {6.00003, 1.64722140125},
     {8.00801, 1.73314962},     {10.00423, 1.77774368},    {12.00499, 1.7941673025},
     {13.99973, 1.8211574975},  {16.99844, 1.8492625775},  {19.99196, 1.869211125},
     {24.98451, 1.8937012575},  {30.00111, 1.9118043575},  {29.99774, 1.9118043575},
     {39.99946, 1.9377028225},  {40.00171, 1.93792921875}, {50.00823, 1.98927441875},
     {50.0099, 1.9847093375},   {60.01801, 1.99381479},    {60.01334, 1.98884559625},
     {70.06358, 2.007667335},   {70.05199, 2.00148027375}, {80.0974, 2.0165039075},
     {80.06891, 2.011145485},   {90.11391, 2.02294543625}, {90.06784, 2.0177893825},
     {100.1277, 2.0289137525},  {100.0715, 2.0246957675},  {120.0024, 2.0441105725},
     {120.0003, 2.04426222},    {140.0057, 2.05697002},    {140.0057, 2.0591104175},
     {160.0078, 2.07108057625}, {160.0096, 2.07156906},    {180.0094, 2.08515767875},
     {180.0028, 2.084146945},   {200.0146, 2.10379276875}, {200.0078, 2.10196756125},
     {220.0133, 2.1210172725},  {220.0051, 2.12265498625}, {240.0215, 2.15010547625},
     {240.0136, 2.1493071775},  {260.0099, 2.1811145375},  {260.0098, 2.1824932275},
     {280.0201, 2.204527525},   {280.0309, 2.20698494625}, {300.0334, 2.239098315},
     {300.0271, 2.2457692275}};

void run_and_print(
    model::ModelInput model,
    common::physical_optimization::OptimizationList optimizationList) {
    runner::Runner runner(model, std::move(optimizationList));
    runner.initializeExperimentalValues(
        experimental_values,
        magnetic_susceptibility::chiT_in_cm_cubed_kelvin_per_mol,
        1);

    runner.minimizeResidualError(std::make_shared<nonlinear_solver::optimNMAdapter>());

    auto theor_final = runner.getMagneticSusceptibilityController().getTheoreticalValues();

    std::cout << std::endl;
    for (auto [temperature, value] : theor_final) {
        std::cout << temperature << " " << value << std::endl;
    }
}

void nickel_crystall() {
    std::cout << "------------------" << std::endl << "Nickel crystal:" << std::endl;

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort().EliminatePositiveProjections().Symmetrize(
        group::Group::S2,
        {{3, 4, 5, 0, 1, 2}});

    // Magnetic motif:
    //    R3b - R3a - Ni
    //    |           |
    //    Ni  - R3a - R3b
    // Spin(Ni) = 1, Spin(R3a) = Spin(R3b) = 1/2

    const std::vector<spin_algebra::Multiplicity> mults = {3, 2, 2, 3, 2, 2};

    // First model.
    std::cout << "\nFirst model\n" << std::endl;

    model::ModelInput first_model(mults);

    const double J_Ni_R3a_value = 40.6;  // K
    auto J_Ni_R3a = first_model.modifySymbolicWorker().addSymbol("J_Ni_R3a", J_Ni_R3a_value);
    first_model.modifySymbolicWorker()
        .assignSymbolToIsotropicExchange(J_Ni_R3a, 0, 1)
        .assignSymbolToIsotropicExchange(J_Ni_R3a, 3, 4);

    const double J_Ni_R3b_value = -333;  // K
    auto J_Ni_R3b = first_model.modifySymbolicWorker().addSymbol("J_Ni_R3b", J_Ni_R3b_value);
    first_model.modifySymbolicWorker()
        .assignSymbolToIsotropicExchange(J_Ni_R3b, 5, 0)
        .assignSymbolToIsotropicExchange(J_Ni_R3b, 2, 3);

    const double J_DR3_value = 4.30;  // K
    auto J_DR3 = first_model.modifySymbolicWorker().addSymbol("J_DR3", J_DR3_value);
    first_model.modifySymbolicWorker()
        .assignSymbolToIsotropicExchange(J_DR3, 1, 2)
        .assignSymbolToIsotropicExchange(J_DR3, 4, 5);

    // R1 and R2 have equal _fixed_ g factor = 2.002
    const double g_R3s_value = 2.002;
    auto g_R3s = first_model.modifySymbolicWorker().addSymbol("g_R3s", g_R3s_value, false);
    first_model.modifySymbolicWorker()
        .assignSymbolToGFactor(g_R3s, 1)
        .assignSymbolToGFactor(g_R3s, 2)
        .assignSymbolToGFactor(g_R3s, 4)
        .assignSymbolToGFactor(g_R3s, 5);
    // Nickels have different free g factor.
    const double g_Ni_value = 2.30;
    auto g_Ni = first_model.modifySymbolicWorker().addSymbol("g_Ni", g_Ni_value);
    first_model.modifySymbolicWorker().assignSymbolToGFactor(g_Ni, 0).assignSymbolToGFactor(
        g_Ni,
        3);

    run_and_print(first_model, optimizationList);

    // J_DR3: 3.56647
    // J_Ni_R3a: 296.944
    // J_Ni_R3b: -322.767
    // g_Ni: 2.01582
    // RSS = 2.92283

    // Second model.
    std::cout << "\nSecond model\n" << std::endl;

    auto second_model = first_model;

    // R2 -- R2 will be replaced by Weiss constant:
    const double Theta_value = 0.0;  // K
    auto Theta = second_model.modifySymbolicWorker().addSymbol("Theta", Theta_value);
    second_model.modifySymbolicWorker().assignSymbolToTheta(Theta);

    run_and_print(second_model, optimizationList);

    // J_DR3: 1.12287
    // J_Ni_R3a: 314.38
    // J_Ni_R3b: -337.919
    // Theta: -1.18885
    // g_Ni: 2.03226
    // RSS = 1.3057

    // Third model.

    std::cout << "\nThird model\n" << std::endl;

    auto third_model = first_model;

    const double D_Ni_value = 6.0;  // K
    auto D = third_model.modifySymbolicWorker().addSymbol("D", D_Ni_value);
    third_model.modifySymbolicWorker()
        .assignSymbolToZFSNoAnisotropy(D, 0)
        .assignSymbolToZFSNoAnisotropy(D, 3);

    run_and_print(third_model, optimizationList);

    // D: 6.72225
    // J_DR3: 2.54587
    // J_Ni_R3a: 307.766
    // J_Ni_R3b: -330.667
    // g_Ni: 2.02476
    // RSS = 1.68627
}
