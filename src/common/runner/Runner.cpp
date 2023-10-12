#include "Runner.h"

#include <cassert>
#include <utility>

#include "src/eigendecompositor/ExactEigendecompositor.h"
#include "src/eigendecompositor/ExplicitQuantitiesEigendecompositor.h"
#include "src/eigendecompositor/ImplicitSSquareEigendecompositor.h"
#include "src/eigendecompositor/OneSymbolInHamiltonianEigendecompositor.h"
#include "src/entities/magnetic_susceptibility/worker/CurieWeissWorker.h"
#include "src/entities/magnetic_susceptibility/worker/GSzSquaredWorker.h"
#include "src/entities/magnetic_susceptibility/worker/UniqueGOnlySSquaredWorker.h"
#include "src/space/optimization/OptimizedSpaceConstructor.h"

namespace runner {

Runner::Runner(model::ModelInput model) :
    Runner(
        std::move(model),
        common::physical_optimization::OptimizationList(),
        quantum::linear_algebra::FactoriesList()) {}

Runner::Runner(
    model::ModelInput model,
    common::physical_optimization::OptimizationList optimizationList) :
    Runner(
        std::move(model),
        std::move(optimizationList),
        quantum::linear_algebra::FactoriesList()) {}

Runner::Runner(
    model::ModelInput model,
    quantum::linear_algebra::FactoriesList dataStructuresFactories) :
    Runner(
        std::move(model),
        common::physical_optimization::OptimizationList(),
        std::move(dataStructuresFactories)) {}

Runner::Runner(
    model::ModelInput model,
    common::physical_optimization::OptimizationList optimizationList,
    quantum::linear_algebra::FactoriesList dataStructuresFactories) :
    consistentModelOptimizationList_(std::move(model), std::move(optimizationList)),
    dataStructuresFactories_(std::move(dataStructuresFactories)),
    space_(space::optimization::OptimizedSpaceConstructor::construct(
        consistentModelOptimizationList_,
        dataStructuresFactories)) {
    // TODO: move it from Runner
    std::unique_ptr<eigendecompositor::AbstractEigendecompositor> eigendecompositor =
        std::make_unique<eigendecompositor::ExactEigendecompositor>(
            getIndexConverter(),
            getDataStructuresFactories());

    // todo: we need only J and D if there is no field
    size_t number_of_changeable_J =
        getModel().getSymbolicWorker().getChangeableNames(model::symbols::J).size();
    size_t number_of_changeable_D =
        getModel().getSymbolicWorker().getChangeableNames(model::symbols::D).size();
    if (number_of_changeable_J == 1 && number_of_changeable_D == 0
        || number_of_changeable_J == 0 && number_of_changeable_D == 1) {
        model::symbols::SymbolName symbol_name;
        if (number_of_changeable_J == 1) {
            symbol_name = getModel().getSymbolicWorker().getChangeableNames(model::symbols::J)[0];
        } else if (number_of_changeable_D == 1) {
            symbol_name = getModel().getSymbolicWorker().getChangeableNames(model::symbols::D)[0];
        }
        auto getter = [this, symbol_name]() {
            return getModel().getSymbolicWorker().getValueOfName(symbol_name);
        };
        eigendecompositor =
            std::make_unique<eigendecompositor::OneSymbolInHamiltonianEigendecompositor>(
                std::move(eigendecompositor),
                getter);
    }

    if (consistentModelOptimizationList_.getOptimizationList().isSSquaredTransformed()) {
        eigendecompositor = std::make_unique<eigendecompositor::ImplicitSSquareEigendecompositor>(
            std::move(eigendecompositor),
            dataStructuresFactories_);
    }

    eigendecompositor = std::make_unique<eigendecompositor::ExplicitQuantitiesEigendecompositor>(
        std::move(eigendecompositor),
        getIndexConverter(),
        getDataStructuresFactories());

    eigendecompositor_ = std::move(eigendecompositor);
}

const space::Space& runner::Runner::getSpace() const {
    return space_;
}

void Runner::BuildSpectra() {
    eigendecompositor_->BuildSpectra(
        consistentModelOptimizationList_.getOperatorsForExplicitConstruction(),
        consistentModelOptimizationList_.getDerivativeOperatorsForExplicitConstruction(),
        getSpace());
}

std::optional<std::reference_wrapper<const Matrix>>
Runner::getMatrix(common::QuantityEnum quantity_enum) const {
    return eigendecompositor_->getMatrix(quantity_enum);
}

const Spectrum& Runner::getSpectrum(common::QuantityEnum quantity_enum) const {
    return eigendecompositor_->getSpectrum(quantity_enum)->get();
}

std::optional<std::shared_ptr<const model::operators::Operator>>
Runner::getOperator(common::QuantityEnum quantity_enum) const {
    return getModel().getOperator(quantity_enum);
}

const lexicographic::IndexConverter& Runner::getIndexConverter() const {
    return getModel().getIndexConverter();
}

std::optional<std::shared_ptr<const model::operators::Operator>> Runner::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return getModel().getOperatorDerivative(quantity_enum, symbol);
}

const Spectrum& Runner::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return eigendecompositor_->getSpectrumDerivative(quantity_enum, symbol)->get();
}

std::optional<std::reference_wrapper<const Matrix>> Runner::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    const model::symbols::SymbolName& symbol) const {
    return eigendecompositor_->getMatrixDerivative(quantity_enum, symbol);
}

void Runner::BuildMuSquaredWorker() {
    auto energy_vector = dataStructuresFactories_.createVector();
    auto degeneracy_vector = dataStructuresFactories_.createVector();

    for (const auto& subspectrum : getSpectrum(common::Energy).blocks) {
        energy_vector->concatenate_with(subspectrum.raw_data);
        degeneracy_vector->add_identical_values(
            subspectrum.raw_data->size(),
            subspectrum.properties.degeneracy);
    }
    energy_vector->subtract_minimum();

    std::unique_ptr<magnetic_susceptibility::worker::AbstractWorker> magnetic_susceptibility_worker;

    if (getSymbolicWorker().isAllGFactorsEqual() && !getSymbolicWorker().isZFSInitialized()) {
        // and there is no field
        // TODO: avoid using of .at(). Change isAllGFactorsEqual signature?
        double g_factor = getModel().getNumericalWorker().getGFactorParameters()->at(0);
        auto s_squared_vector = dataStructuresFactories_.createVector();

        for (const auto& subspectrum : getSpectrum(common::S_total_squared).blocks) {
            s_squared_vector->concatenate_with(subspectrum.raw_data);
        }

        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::UniqueGOnlySSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(s_squared_vector),
                g_factor);
    } else {
        auto g_sz_squared_vector = dataStructuresFactories_.createVector();
        // TODO: check if g_sz_squared has been initialized
        for (const auto& subspectrum : getSpectrum(common::gSz_total_squared).blocks) {
            g_sz_squared_vector->concatenate_with(subspectrum.raw_data);
        }

        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::GSzSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(g_sz_squared_vector));
    }

    if (getSymbolicWorker().isThetaInitialized()) {
        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::CurieWeissWorker>(
                std::move(magnetic_susceptibility_worker),
                getModel().getNumericalWorker().getThetaParameter());
    }

    magnetic_susceptibility_controller_ = magnetic_susceptibility::MagneticSusceptibilityController(
        std::move(magnetic_susceptibility_worker));

    if (experimental_values_worker_.has_value()) {
        magnetic_susceptibility_controller_.value().initializeExperimentalValues(
            experimental_values_worker_.value());
    }
}

void Runner::initializeExperimentalValues(
    const std::vector<magnetic_susceptibility::ValueAtTemperature>& experimental_data,
    magnetic_susceptibility::ExperimentalValuesEnum experimental_quantity_type,
    double number_of_centers_ratio) {
    if (experimental_values_worker_.has_value()) {
        throw std::invalid_argument("Experimental values have been already initialized");
    }

    experimental_values_worker_ =
        std::make_shared<magnetic_susceptibility::ExperimentalValuesWorker>(
            experimental_data,
            experimental_quantity_type,
            number_of_centers_ratio);

    if (magnetic_susceptibility_controller_.has_value()) {
        magnetic_susceptibility_controller_.value().initializeExperimentalValues(
            experimental_values_worker_.value());
    }
}

std::map<model::symbols::SymbolName, double> Runner::calculateTotalDerivatives() {
    // TODO: it is awful. Fix it somehow, this code should be moved from Runner.

    std::map<model::symbols::SymbolName, double> answer;

    for (const auto& changeable_symbol :
         getSymbolicWorker().getChangeableNames(model::symbols::J)) {
        auto derivative_vector = dataStructuresFactories_.createVector();
        for (const auto& subspectrum :
             getSpectrumDerivative(common::Energy, changeable_symbol).blocks) {
            derivative_vector->concatenate_with(subspectrum.raw_data);
        }
        auto derivative_map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        derivative_map[common::Energy] = std::move(derivative_vector);
        double value = magnetic_susceptibility_controller_.value().calculateTotalDerivative(
            model::symbols::J,
            std::move(derivative_map));
        answer[changeable_symbol] = value;
        //        std::cout << "dRSS/d" << changeable_symbol.get_name() << " = " << value << std::endl;
    }

    for (const auto& changeable_symbol :
         getSymbolicWorker().getChangeableNames(model::symbols::D)) {
        auto derivative_vector = dataStructuresFactories_.createVector();
        for (const auto& subspectrum :
             getSpectrumDerivative(common::Energy, changeable_symbol).blocks) {
            derivative_vector->concatenate_with(subspectrum.raw_data);
        }
        auto derivative_map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        derivative_map[common::Energy] = std::move(derivative_vector);
        double value = magnetic_susceptibility_controller_.value().calculateTotalDerivative(
            model::symbols::D,
            std::move(derivative_map));
        answer[changeable_symbol] = value;
        //        std::cout << "dRSS/d" << changeable_symbol.get_name() << " = " << value << std::endl;
    }

    for (const auto& changeable_symbol :
         getSymbolicWorker().getChangeableNames(model::symbols::g_factor)) {
        auto map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        if (consistentModelOptimizationList_.getModel().is_g_sz_squared_derivatives_initialized()) {
            auto derivative_vector = dataStructuresFactories_.createVector();
            for (const auto& subspectrum :
                 getSpectrumDerivative(common::gSz_total_squared, changeable_symbol).blocks) {
                derivative_vector->concatenate_with(subspectrum.raw_data);
            }
            map[common::gSz_total_squared] = std::move(derivative_vector);
        }
        double value = magnetic_susceptibility_controller_.value().calculateTotalDerivative(
            model::symbols::g_factor,
            std::move(map));
        answer[changeable_symbol] = value;
        //        std::cout << "dRSS/d" << changeable_symbol.get_name() << " = " << value << std::endl;
    }

    // Theta calculation:
    if (!getSymbolicWorker().getChangeableNames(model::symbols::Theta).empty()) {
        model::symbols::SymbolName Theta_name =
            getSymbolicWorker().getChangeableNames(model::symbols::Theta)[0];
        auto empty_map = std::map<
            common::QuantityEnum,
            std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>();
        double value = magnetic_susceptibility_controller_.value().calculateTotalDerivative(
            model::symbols::Theta,
            std::move(empty_map));
        answer[Theta_name] = value;
        //        std::cout << "dRSS/d" << Theta_name.get_name() << " = " << value << std::endl;
    }

    return answer;
}

void Runner::minimizeResidualError(
    std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver> solver) {
    std::vector<model::symbols::SymbolName> changeable_names =
        getSymbolicWorker().getChangeableNames();
    std::vector<double> changeable_values;
    changeable_values.reserve(changeable_names.size());
    for (const model::symbols::SymbolName& name : changeable_names) {
        changeable_values.push_back(getSymbolicWorker().getValueOfName(name));
    }

    if (solver->doesGradientsRequired()) {
        initializeDerivatives();
    }

    // This function should calculate residual error and derivatives at a point of changeable_values
    std::function<double(const std::vector<double>&, std::vector<double>&, bool)> oneStepFunction =
        [this, capture0 = std::cref(changeable_names)](
            const std::vector<double>& changeable_values,
            std::vector<double>& gradient,
            bool isGradientRequired) {
            return stepOfRegression(capture0, changeable_values, gradient, isGradientRequired);
        };

    solver->optimize(oneStepFunction, changeable_values);

    for (size_t i = 0; i < changeable_names.size(); ++i) {
        std::cout << changeable_names[i].get_name() << ": " << changeable_values[i] << std::endl;
    }
    std::cout << "RSS = " << magnetic_susceptibility_controller_.value().calculateResidualError()
              << std::endl;
}

double Runner::stepOfRegression(
    const std::vector<model::symbols::SymbolName>& changeable_names,
    const std::vector<double>& changeable_values,
    std::vector<double>& gradient,
    bool isGradientRequired) {
    // At first, update actual values in SymbolicWorker:
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        consistentModelOptimizationList_.setNewValueToChangeableSymbol(
            changeable_names[i],
            changeable_values[i]);
    }

    //    for (size_t i = 0; i < changeable_names.size(); ++i) {
    //        std::cout << changeable_names[i].get_name() << " = " << changeable_values[i] << std::endl;
    //    }

    // Do some calculation stuff...
    BuildSpectra();

    BuildMuSquaredWorker();

    // Calculate residual error and write it to external variable:
    double residual_error = getMagneticSusceptibilityController().calculateResidualError();

    std::cout << "RSS = " << residual_error << std::endl << std::endl;

    if (isGradientRequired) {
        // Calculate derivatives...
        std::map<model::symbols::SymbolName, double> map_gradient = calculateTotalDerivatives();
        for (size_t i = 0; i < changeable_names.size(); ++i) {
            // ...and write it to external variable:
            gradient[i] = map_gradient[changeable_names[i]];
        }
    }

    return residual_error;
}

void Runner::initializeDerivatives() {
    consistentModelOptimizationList_.InitializeDerivatives();
}

const magnetic_susceptibility::MagneticSusceptibilityController&
Runner::getMagneticSusceptibilityController() const {
    return magnetic_susceptibility_controller_.value();
}

const model::symbols::SymbolicWorker& Runner::getSymbolicWorker() const {
    return getModel().getSymbolicWorker();
}

const model::Model& Runner::getModel() const {
    return consistentModelOptimizationList_.getModel();
}

const common::physical_optimization::OptimizationList& Runner::getOptimizationList() const {
    return consistentModelOptimizationList_.getOptimizationList();
}

quantum::linear_algebra::FactoriesList Runner::getDataStructuresFactories() const {
    return dataStructuresFactories_;
}
}  // namespace runner