#include "src/common/runner/Runner.h"
#include "src/nonlinear_solver/optimNM/optimNMAdapter.h"

// Calculation of chiT of manganese cocrystal from https://doi.org/10.1039/D0DT03184D
void manganese_cocrystall() {
    std::cout << "Manganese cocrystal:" << std::endl;
    // Values: cm^3K/mol
    std::vector<magnetic_susceptibility::ValueAtTemperature> values = {
        {1.99976, 1.49992},  {1.99921, 1.49874},  {2.49982, 1.74596},  {2.99904, 1.95376},
        {3.49923, 2.12467},  {4.01445, 2.26855},  {4.49826, 2.39634},  {4.99923, 2.51065},
        {6.00094, 2.69048},  {8.00233, 2.91901},  {10.00243, 3.0579},  {12.00708, 3.13307},
        {14.00121, 3.21003}, {16.99602, 3.29},    {19.99461, 3.34786}, {24.98512, 3.41624},
        {29.99732, 3.46265}, {29.99015, 3.46125}, {39.99637, 3.52593}, {39.99553, 3.52467},
        {50.00529, 3.60335}, {49.99944, 3.60836}, {60.04366, 3.6842},  {70.07054, 3.67018},
        {70.02648, 3.65688}, {80.06694, 3.68126}, {80.06104, 3.67168}, {90.06602, 3.70037},
        {100.0996, 3.74245}, {100.0575, 3.73534}, {120.0024, 3.83788}, {120.0045, 3.83747},
        {140.0002, 3.96887}, {140.0002, 3.96845}, {160.0096, 4.12135}, {160.0113, 4.12164},
        {180.0094, 4.29182}, {180.0116, 4.29251}, {200.0132, 4.47581}, {200.0105, 4.47584},
        {220.0149, 4.66605}, {220.0116, 4.66656}, {240.0039, 4.86564}, {240.0156, 4.86623},
        {260.0213, 5.08483}, {260.0122, 5.08365}, {280.0255, 5.27514}, {280.0201, 5.27526},
        {300.0085, 5.4826},  {299.9866, 5.48317}};
    std::vector<int> mults = {6, 2, 2, 6, 2, 2};
    model::ModelInput model(mults);

    // Magnetic motif:
    // -- R2 - R1 - Mn
    //    |         |
    //    Mn - R1 - R2 --
    // Spin(Mn) = 5/2, Spin(R1) = Spin(R2) = 1/2

    // It is cocrystal of two Mn complexes with proportion 3:2.
    // So we choose guess values as weighted average of these two complexes values.

    // (3*(-103)+2*(-43))/5=-79
    double J_Mn_R1_initial = -79.0;
    auto J_Mn_R1 = model.modifySymbolicWorker().addSymbol("J_Mn_R1", J_Mn_R1_initial);
    model.modifySymbolicWorker()
        .assignSymbolToIsotropicExchange(J_Mn_R1, 0, 1)
        .assignSymbolToIsotropicExchange(J_Mn_R1, 3, 4);
    // (3*(-125)+2*(-114))/5=-120.6
    double J_Mn_R2_initial = -120.6;
    auto J_Mn_R2 = model.modifySymbolicWorker().addSymbol("J_Mn_R2", J_Mn_R2_initial);
    model.modifySymbolicWorker()
        .assignSymbolToIsotropicExchange(J_Mn_R2, 5, 0)
        .assignSymbolToIsotropicExchange(J_Mn_R2, 2, 3);
    // (3*9.9+2*28.1)/5=17.18
    double J_dirad_initial = 17.18;
    auto J_dirad = model.modifySymbolicWorker().addSymbol("J_dirad", J_dirad_initial);
    model.modifySymbolicWorker()
        .assignSymbolToIsotropicExchange(J_dirad, 1, 2)
        .assignSymbolToIsotropicExchange(J_dirad, 4, 5);
    // All g factors are equal.
    double g_initial = 2.002;
    auto g = model.modifySymbolicWorker().addSymbol("g", g_initial);
    for (size_t i = 0; i < mults.size(); ++i) {
        model.modifySymbolicWorker().assignSymbolToGFactor(g, i);
    }
    // R2 -- R2 will be replaced by Weiss constant:
    // (3*(-11.7)+2*(-11.5))/5=-11.62
    double Theta_initial = -11.62;
    auto Theta = model.modifySymbolicWorker().addSymbol("Theta", Theta_initial);
    model.modifySymbolicWorker().assignSymbolToTheta(Theta);

    common::physical_optimization::OptimizationList optimizationList;
    optimizationList.TzSort().EliminatePositiveProjections().Symmetrize(
        group::Group::S2,
        {{3, 4, 5, 0, 1, 2}});

    runner::Runner runner(model, optimizationList);
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
    // J_Mn_R1: -111.46
    // J_Mn_R2: -89.6323
    // J_dirad: 16.382
    // Theta: -5.99865
    // g: 2.00861
    // RSS = 0.682388

    // It should be mentioned, that J_Mn_R1 and J_Mn_R2 swapped during optimization
    // (it's okay due to symmetry of Hamiltonian). Thus, the final answer:
    //
    // J_Mn_R1 = -89.6
    // J_Mn_R2 = -111
    // J_dirad =  16.4
    // Theta   = -6.00
    // g:      =  2.009
    // RSS     =  0.682388
}
