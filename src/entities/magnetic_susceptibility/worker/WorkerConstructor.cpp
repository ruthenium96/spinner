#include "WorkerConstructor.h"

#include "src/common/Logger.h"
#include "src/common/Quantity.h"
#include "src/entities/magnetic_susceptibility/worker/CurieWeissWorker.h"
#include "src/entities/magnetic_susceptibility/worker/GSzSquaredWorker.h"
#include "src/entities/magnetic_susceptibility/worker/UniqueGOnlySSquaredWorker.h"

namespace magnetic_susceptibility::worker {

std::unique_ptr<AbstractWorker> WorkerConstructor::construct(
    const model::Model& model,
    const std::unique_ptr<eigendecompositor::AbstractEigendecompositor>& eigendecompositor,
    const quantum::linear_algebra::FactoriesList& factories) {

    common::Logger::detailed_msg("magnetic_susceptibility::Worker information:");
    
    auto energy_vector = factories.createVector();
    auto degeneracy_vector = factories.createVector();

    for (const auto& subspectrum : eigendecompositor->getSpectrum(common::Energy).value().get().blocks) {
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
             eigendecompositor->getSpectrum(common::S_total_squared).value().get().blocks) {
            s_squared_vector->concatenate_with(subspectrum.raw_data);
        }

        common::Logger::detailed("UniqueGOnlySSquaredWorker will be used.");
        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::UniqueGOnlySSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(s_squared_vector),
                g_factor);
    } else if (model.is_g_sz_squared_initialized()) {
        auto g_sz_squared_vector = factories.createVector();
        for (const auto& subspectrum :
             eigendecompositor->getSpectrum(common::gSz_total_squared).value().get().blocks) {
            g_sz_squared_vector->concatenate_with(subspectrum.raw_data);
        }

        common::Logger::detailed("GSzSquaredWorker will be used.");
        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::GSzSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(g_sz_squared_vector));
    } else {
        throw std::invalid_argument("Cannot construct magnetic_susceptibility::worker");
    }

    if (symbolic_worker.isThetaInitialized()) {
        common::Logger::detailed("CurieWeissWorker will be used.");
        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::CurieWeissWorker>(
                std::move(magnetic_susceptibility_worker),
                model.getNumericalWorker().getThetaParameter());
    }

    common::Logger::separate(0, common::detailed);

    return magnetic_susceptibility_worker;

}
}  // namespace magnetic_susceptibility::worker