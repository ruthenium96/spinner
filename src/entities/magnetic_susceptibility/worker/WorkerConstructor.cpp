#include "WorkerConstructor.h"
#include <stdexcept>

#include "src/common/Logger.h"
#include "src/common/Quantity.h"
#include "src/entities/magnetic_susceptibility/worker/CurieWeissWorker.h"
#include "src/entities/magnetic_susceptibility/worker/DifferentGWorker.h"
#include "src/entities/magnetic_susceptibility/worker/UniqueGWorker.h"

namespace magnetic_susceptibility::worker {

std::unique_ptr<AbstractWorker> WorkerConstructor::construct(
    const runner::ConsistentModelOptimizationList& consistentModelOptimizationList,
    std::shared_ptr<const eigendecompositor::FlattenedSpectra> flattenedSpectra,
    const quantum::linear_algebra::FactoriesList& factories) {

    common::Logger::detailed_msg("magnetic_susceptibility::Worker information:");
    
    std::unique_ptr<magnetic_susceptibility::worker::AbstractWorker> magnetic_susceptibility_worker;

    const auto& symbolic_worker = consistentModelOptimizationList.getModel().getSymbolicWorker();

    if (consistentModelOptimizationList.isImplicitMSquarePossible() 
        || consistentModelOptimizationList.isImplicitSSquarePossible()
        || consistentModelOptimizationList.isExplicitMSquarePossible()) {
        // and there is no field
        // TODO: avoid using of .at(). Change isAllGFactorsEqual signature?
        double g_factor = consistentModelOptimizationList.getModel().getNumericalWorker().getGFactorParameters()->at(0);
        common::QuantityEnum quantity_enum_for_averaging;
        double quantity_factor;
        if (flattenedSpectra->getFlattenSpectrum(common::S_total_squared).has_value()) {
            quantity_enum_for_averaging = common::S_total_squared;
            quantity_factor = 1.0;
            common::Logger::detailed("UniqueGWorker will be used with S_total_squared.");
        } else if (flattenedSpectra->getFlattenSpectrum(common::M_total_squared).has_value()) {
            quantity_enum_for_averaging = common::M_total_squared;
            quantity_factor = 3.0;
            common::Logger::detailed("UniqueGWorker will be used with M_total_squared.");
        } else {
            throw std::invalid_argument("UniqueGWorker cannot find Spectrum of S_total_squared or M_total_squared");
        }
        
        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::UniqueGWorker>(
                flattenedSpectra,
                g_factor, 
                quantity_enum_for_averaging, 
                quantity_factor);
    } else if (consistentModelOptimizationList.isGSquaredT00Possible() 
        || consistentModelOptimizationList.isGSzSquaredPossible()) {
        common::QuantityEnum quantity_enum_for_averaging;
        double quantity_factor;
        if (flattenedSpectra->getFlattenSpectrum(common::g_squared_T00).has_value()) {
            quantity_enum_for_averaging = common::g_squared_T00;
            quantity_factor = 3.0;
            common::Logger::detailed("GSzSquaredWorker will be used with g_squared_T00.");
        } else if (flattenedSpectra->getFlattenSpectrum(common::gSz_total_squared).has_value()) {
            quantity_enum_for_averaging = common::gSz_total_squared;
            quantity_factor = 3.0;
            common::Logger::detailed("GSzSquaredWorker will be used with gSz_total_squared.");
        } else {
            throw std::invalid_argument("GSzSquaredWorker cannot find Spectrum of g_squared_T00 or gSz_total_squared");
        }

        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::DifferentGWorker>(
                flattenedSpectra,
                quantity_enum_for_averaging,
                quantity_factor);
    } else {
        throw std::invalid_argument("Cannot construct magnetic_susceptibility::worker");
    }

    if (symbolic_worker.isThetaInitialized()) {
        common::Logger::detailed("CurieWeissWorker will be used.");
        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::CurieWeissWorker>(
                std::move(magnetic_susceptibility_worker),
                consistentModelOptimizationList.getModel().getNumericalWorker().getThetaParameter());
    }

    common::Logger::separate(0, common::detailed);

    return magnetic_susceptibility_worker;

}
}  // namespace magnetic_susceptibility::worker