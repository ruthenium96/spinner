#include "WorkerConstructor.h"

#include "src/entities/magnetic_susceptibility/worker/CurieWeissWorker.h"
#include "src/entities/magnetic_susceptibility/worker/GSzSquaredWorker.h"
#include "src/entities/magnetic_susceptibility/worker/UniqueGOnlySSquaredWorker.h"

namespace magnetic_susceptibility::worker {

std::unique_ptr<AbstractWorker> WorkerConstructor::construct(
    const model::Model& model,
    const std::unique_ptr<eigendecompositor::AbstractEigendecompositor>& eigendecompositor,
    const quantum::linear_algebra::FactoriesList& factories) {
    auto energy_vector = factories.createVector();
    auto degeneracy_vector = factories.createVector();

    for (const auto& subspectrum : eigendecompositor->getSpectrum(common::Energy)->get().blocks) {
        energy_vector->concatenate_with(subspectrum.raw_data);
        degeneracy_vector->add_identical_values(
            subspectrum.raw_data->size(),
            subspectrum.properties.degeneracy);
    }
    energy_vector->subtract_minimum();

    std::unique_ptr<magnetic_susceptibility::worker::AbstractWorker> magnetic_susceptibility_worker;

    const auto& symbolic_worker = model.getSymbolicWorker();

    if (symbolic_worker.isAllGFactorsEqual() && !symbolic_worker.isZFSInitialized()) {
        // and there is no field
        // TODO: avoid using of .at(). Change isAllGFactorsEqual signature?
        double g_factor = model.getNumericalWorker().getGFactorParameters()->at(0);
        auto s_squared_vector = factories.createVector();

        for (const auto& subspectrum :
             eigendecompositor->getSpectrum(common::S_total_squared)->get().blocks) {
            s_squared_vector->concatenate_with(subspectrum.raw_data);
        }

        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::UniqueGOnlySSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(s_squared_vector),
                g_factor);
    } else {
        auto g_sz_squared_vector = factories.createVector();
        // TODO: check if g_sz_squared has been initialized
        for (const auto& subspectrum :
             eigendecompositor->getSpectrum(common::gSz_total_squared)->get().blocks) {
            g_sz_squared_vector->concatenate_with(subspectrum.raw_data);
        }

        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::GSzSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(g_sz_squared_vector));
    }

    if (symbolic_worker.isThetaInitialized()) {
        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::CurieWeissWorker>(
                std::move(magnetic_susceptibility_worker),
                model.getNumericalWorker().getThetaParameter());
    }

    return magnetic_susceptibility_worker;

}
}  // namespace magnetic_susceptibility::worker