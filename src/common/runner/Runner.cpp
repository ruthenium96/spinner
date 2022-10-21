#include "Runner.h"

#include <cassert>
#include <utility>

#include "src/model/operators/ConstantTerm.h"
#include "src/model/operators/ScalarProductTerm.h"
#include "src/space/optimization/OptimizedSpaceConstructor.h"

namespace runner {

Runner::Runner(model::Model model) :
    Runner(
        std::move(model),
        common::physical_optimization::OptimizationList(),
        quantum::linear_algebra::AbstractFactory::defaultFactory()) {}

Runner::Runner(
    model::Model model,
    common::physical_optimization::OptimizationList optimizationList) :
    Runner(
        std::move(model),
        std::move(optimizationList),
        quantum::linear_algebra::AbstractFactory::defaultFactory()) {}

Runner::Runner(
    model::Model model,
    std::shared_ptr<quantum::linear_algebra::AbstractFactory> algebraDataFactory) :
    Runner(
        std::move(model),
        common::physical_optimization::OptimizationList(),
        std::move(algebraDataFactory)) {}

Runner::Runner(
    model::Model model,
    common::physical_optimization::OptimizationList optimizationList,
    std::shared_ptr<quantum::linear_algebra::AbstractFactory> algebraDataFactory) :
    consistentModelOptimizationList_(std::move(model), std::move(optimizationList)),
    algebraDataFactory_(std::move(algebraDataFactory)),
    space_(space::optimization::OptimizedSpaceConstructor::construct(
        consistentModelOptimizationList_)) {
    if (getModel().is_s_squared_initialized()) {
        s_squared = common::Quantity();
    }
    if (getModel().is_isotropic_exchange_derivatives_initialized()) {
        for (const auto& symbol :
             getSymbols().getChangeableNames(model::symbols::SymbolTypeEnum::J)) {
            derivative_of_energy_wrt_exchange_parameters[symbol] = common::Quantity();
        }
    }

    if (getSymbols().isIsotropicExchangeInitialized()) {
        getModel().InitializeIsotropicExchange();
    }

    //    if (!symbols_.isGFactorInitialized()) {
    //        throw std::length_error("g factor parameters have not been initialized");
    //    }
}

const space::Space& runner::Runner::getSpace() const {
    return space_;
}

void Runner::BuildMatrices() {
    if (!getOperator(common::Energy).empty()) {
        energy.matrix_ = Matrix(
            getSpace(),
            getOperator(common::Energy),
            getIndexConverter(),
            getAlgebraDataFactory());
    }
    if (s_squared.has_value()) {
        s_squared->matrix_ = Matrix(
            getSpace(),
            getOperator(common::S_total_squared),
            getIndexConverter(),
            getAlgebraDataFactory());
    }
    for (auto& [symbol_name, derivative] : derivative_of_energy_wrt_exchange_parameters) {
        // TODO: fix it!
        derivative.matrix_ = Matrix(
            getSpace(),
            getOperatorDerivative(common::Energy, model::symbols::J, symbol_name),
            getIndexConverter(),
            getAlgebraDataFactory());
    }

    matrix_history_.matrices_was_built = true;
}

void Runner::BuildSpectra() {
    size_t number_of_blocks = getSpace().getBlocks().size();

    if (!getOperator(common::Energy).empty() || !energy.spectrum_.blocks.empty()) {
        energy.spectrum_.blocks.clear();
        energy.spectrum_.blocks.reserve(number_of_blocks);
    }
    if (s_squared.has_value()) {
        s_squared->spectrum_.blocks.clear();
        s_squared->spectrum_.blocks.reserve(number_of_blocks);
    }
    for (auto& [_, derivative] : derivative_of_energy_wrt_exchange_parameters) {
        derivative.spectrum_.blocks.clear();
        derivative.spectrum_.blocks.reserve(number_of_blocks);
    }

    if (matrix_history_.matrices_was_built) {
        BuildSpectraUsingMatrices(number_of_blocks);
    } else {
        BuildSpectraWithoutMatrices(number_of_blocks);
    }
}

void Runner::BuildSpectraUsingMatrices(size_t number_of_blocks) {
    for (size_t block = 0; block < number_of_blocks; ++block) {
        auto [subspectrum_energy, unitary_transformation_matrix] =
            Subspectrum::energy(energy.matrix_.blocks[block]);
        energy.spectrum_.blocks.emplace_back(std::move(subspectrum_energy));

        if (s_squared.has_value()) {
            s_squared->spectrum_.blocks.emplace_back(Subspectrum::non_energy(
                s_squared->matrix_.blocks[block],
                unitary_transformation_matrix));
        }

        for (auto& [_, derivative] : derivative_of_energy_wrt_exchange_parameters) {
            derivative.spectrum_.blocks.emplace_back(Subspectrum::non_energy(
                derivative.matrix_.blocks[block],
                unitary_transformation_matrix));
        }
    }
}

void Runner::BuildSpectraWithoutMatrices(size_t number_of_blocks) {
    for (size_t block = 0; block < number_of_blocks; ++block) {
        auto unitary_transformation_matrix = algebraDataFactory_->createMatrix();
        {
            auto hamiltonian_submatrix = Submatrix(
                getSpace().getBlocks()[block],
                getOperator(common::Energy),
                getIndexConverter(),
                getAlgebraDataFactory());
            auto pair = Subspectrum::energy(hamiltonian_submatrix);
            energy.spectrum_.blocks.emplace_back(std::move(pair.first));
            unitary_transformation_matrix = std::move(pair.second);
        }

        if (s_squared.has_value()) {
            auto non_hamiltonian_submatrix = Submatrix(
                getSpace().getBlocks()[block],
                getOperator(common::S_total_squared),
                getIndexConverter(),
                getAlgebraDataFactory());
            s_squared->spectrum_.blocks.emplace_back(
                Subspectrum::non_energy(non_hamiltonian_submatrix, unitary_transformation_matrix));
        }

        for (auto& [symbol_name, derivative] : derivative_of_energy_wrt_exchange_parameters) {
            // TODO: fix it
            auto derivative_submatrix = Submatrix(
                getSpace().getBlocks()[block],
                getOperatorDerivative(common::Energy, model::symbols::J, symbol_name),
                getIndexConverter(),
                getAlgebraDataFactory());
            derivative.spectrum_.blocks.emplace_back(
                Subspectrum::non_energy(derivative_submatrix, unitary_transformation_matrix));
        }
    }
}

const Matrix& Runner::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return energy.matrix_;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return s_squared->matrix_;
    }
    assert(0);
}

const Spectrum& Runner::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return energy.spectrum_;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return s_squared->spectrum_;
    }
    assert(0);
}

const model::operators::Operator& Runner::getOperator(common::QuantityEnum quantity_enum) const {
    return getModel().getOperator(quantity_enum);
}

const lexicographic::IndexConverter& Runner::getIndexConverter() const {
    return getModel().getIndexConverter();
}

const model::operators::Operator& Runner::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    model::symbols::SymbolTypeEnum symbol_type,
    const model::symbols::SymbolName& symbol) const {
    return getModel().getOperatorDerivative(quantity_enum, symbol_type, symbol);
}

const Spectrum& Runner::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    model::symbols::SymbolTypeEnum symbol_type,
    const model::symbols::SymbolName& symbol) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        if (symbol_type == model::symbols::SymbolTypeEnum::J) {
            return derivative_of_energy_wrt_exchange_parameters.at(symbol).spectrum_;
        }
    }
}

const Matrix& Runner::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    model::symbols::SymbolTypeEnum symbol_type,
    const model::symbols::SymbolName& symbol) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        if (symbol_type == model::symbols::SymbolTypeEnum::J) {
            return derivative_of_energy_wrt_exchange_parameters.at(symbol).matrix_;
        }
    }
}

void Runner::BuildMuSquaredWorker() {
    auto energy_vector = algebraDataFactory_->createVector();
    auto degeneracy_vector = algebraDataFactory_->createVector();

    for (const auto& subspectrum : energy.spectrum_.blocks) {
        energy_vector->concatenate_with(subspectrum.raw_data);
        degeneracy_vector->add_identical_values(
            subspectrum.raw_data->size(),
            subspectrum.properties.degeneracy);
    }
    energy_vector->subtract_minimum();

    std::unique_ptr<magnetic_susceptibility::worker::AbstractWorker> magnetic_susceptibility_worker;

    if (getSymbols().isAllGFactorsEqual()) {
        // and there is no field
        // TODO: avoid using of .at(). Change isAllGFactorsEqual signature?
        double g_factor = getSymbols().getGFactorParameters()->at(0);
        auto s_squared_vector = algebraDataFactory_->createVector();

        // TODO: check if s_squared has been initialized
        for (const auto& subspectrum : s_squared->spectrum_.blocks) {
            s_squared_vector->concatenate_with(subspectrum.raw_data);
        }

        magnetic_susceptibility_worker =
            std::make_unique<magnetic_susceptibility::worker::UniqueGOnlySSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(s_squared_vector),
                g_factor);

        if (getSymbols().isThetaInitialized()) {
            magnetic_susceptibility_worker =
                std::make_unique<magnetic_susceptibility::worker::CurieWeissWorker>(
                    std::move(magnetic_susceptibility_worker),
                    getSymbols().getThetaParameter());
        }
    } else {
        throw std::invalid_argument("Different g factors are not supported now.");
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
    // TODO: ony s_squared-based calculation supported
    // TODO: it is awful. Fix it somehow, this code should be moved from Runner.

    std::map<model::symbols::SymbolName, double> answer;

    for (const auto& changeable_symbol : getSymbols().getChangeableNames(model::symbols::J)) {
        auto derivative_vector = algebraDataFactory_->createVector();
        for (const auto& subspectrum :
             getSpectrumDerivative(common::Energy, model::symbols::J, changeable_symbol).blocks) {
            derivative_vector->concatenate_with(subspectrum.raw_data);
        }
        double value = magnetic_susceptibility_controller_.value().calculateTotalDerivative(
            model::symbols::J,
            std::move(derivative_vector));
        answer[changeable_symbol] = value;
        //        std::cout << "dR^2/d" << changeable_symbol.get_name() << " = " << value << std::endl;
    }

    if (!getSymbols().getChangeableNames(model::symbols::g_factor).empty()) {
        model::symbols::SymbolName g_name =
            getSymbols().getChangeableNames(model::symbols::g_factor)[0];
        double value = magnetic_susceptibility_controller_.value().calculateTotalDerivative(
            model::symbols::g_factor);
        answer[g_name] = value;
        //        std::cout << "dR^2/d" << g_name.get_name() << " = " << value << std::endl;
    }

    // Theta calculation:
    if (!getSymbols().getChangeableNames(model::symbols::Theta).empty()) {
        model::symbols::SymbolName Theta_name =
            getSymbols().getChangeableNames(model::symbols::Theta)[0];
        double value = magnetic_susceptibility_controller_.value().calculateTotalDerivative(
            model::symbols::Theta);
        answer[Theta_name] = value;
        //        std::cout << "dR^2/d" << Theta_name.get_name() << " = " << value << std::endl;
    }

    return answer;
}

void Runner::minimizeResidualError(
    std::shared_ptr<nonlinear_solver::AbstractNonlinearSolver> solver) {
    std::vector<model::symbols::SymbolName> changeable_names = getSymbols().getChangeableNames();
    std::vector<double> changeable_values;
    changeable_values.reserve(changeable_names.size());
    for (const model::symbols::SymbolName& name : changeable_names) {
        changeable_values.push_back(getSymbols().getValueOfName(name));
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
        //        std::cout << changeable_names[i].get_name() << ": " << changeable_values[i] << std::endl;
    }
}

double Runner::stepOfRegression(
    const std::vector<model::symbols::SymbolName>& changeable_names,
    const std::vector<double>& changeable_values,
    std::vector<double>& gradient,
    bool isGradientRequired) {
    // At first, update actual values in Symbols:
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        // TODO: mutable use of Model/Symbols. Refactor it.
        getModel().getSymbols().setNewValueToChangeableSymbol(
            changeable_names[i],
            changeable_values[i]);
    }

    // Do some calculation stuff...
    if (matrix_history_.matrices_was_built) {
        // TODO: 1) it does not work now
        //       2) here's the problem with future SSquaredTransformer
        //       or we can use SSquaredTransform() as flag and actually apply it at BuildMatrix?
        BuildMatrices();
    }
    BuildSpectra();

    BuildMuSquaredWorker();

    // Calculate residual error and write it to external variable:
    double residual_error = magnetic_susceptibility_controller_.value().calculateResidualError();

    //    std::cout << "R^2 = " << residual_error << std::endl;

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

const magnetic_susceptibility::MagneticSusceptibilityController&
Runner::getMagneticSusceptibilityController() const {
    return magnetic_susceptibility_controller_.value();
}

const model::symbols::Symbols& Runner::getSymbols() const {
    return getModel().getSymbols();
}

const model::Model& Runner::getModel() const {
    return consistentModelOptimizationList_.getModel();
}

model::Model& Runner::getModel() {
    return consistentModelOptimizationList_.getModel();
}

const common::physical_optimization::OptimizationList& Runner::getOptimizationList() const {
    return consistentModelOptimizationList_.getOptimizationList();
}

std::shared_ptr<quantum::linear_algebra::AbstractFactory> Runner::getAlgebraDataFactory() const {
    return algebraDataFactory_;
}
}  // namespace runner