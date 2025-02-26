#include "WorkerConstructor.h"
#include <stdexcept>

#include "src/common/Logger.h"
#include "src/common/Quantity.h"
#include "src/entities/magnetic_susceptibility/worker/CurieWeissWorker.h"
#include "src/entities/magnetic_susceptibility/worker/GSzSquaredWorker.h"
#include "src/entities/magnetic_susceptibility/worker/UniqueGWorker.h"

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
        auto quantity_vector = factories.createVector();

        if (eigendecompositor->getSpectrum(common::S_total_squared).has_value()) {
            for (const auto& subspectrum :
                eigendecompositor->getSpectrum(common::S_total_squared).value().get().blocks) {
                    quantity_vector->concatenate_with(subspectrum.raw_data);
            }   
            common::Logger::detailed("UniqueGOnlySSquaredWorker will be used with S_total_squared.");
        } else if (eigendecompositor->getSpectrum(common::M_total_squared).has_value()) {
            for (const auto& subspectrum :
                eigendecompositor->getSpectrum(common::M_total_squared).value().get().blocks) {
                    quantity_vector->concatenate_with(subspectrum.raw_data);
            }
            // TODO: create in-place multiply_by function?
            quantity_vector = quantity_vector->multiply_by(3);
            common::Logger::detailed("UniqueGOnlySSquaredWorker will be used with M_total_squared.");
        } else {
            throw std::invalid_argument("UniqueGOnlySSquaredWorker cannot find Spectrum of S_total_squared or M_total_squared");
        }

        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::UniqueGWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(quantity_vector),
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