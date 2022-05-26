#include "Runner.h"

#include <stlbfgs.h>

#include <utility>

#include "src/model/operators/ConstantTerm.h"
#include "src/model/operators/ScalarProductTerm.h"
#include "src/space/optimization/OptimizedSpaceConstructor.h"

namespace runner {

Runner::Runner(model::Model model) :
    Runner(std::move(model), common::physical_optimization::OptimizationList()) {}

Runner::Runner(
    model::Model model,
    common::physical_optimization::OptimizationList optimizationList) :
    consistentModelOptimizationList_(std::move(model), std::move(optimizationList)),
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
        energy.matrix_ = Matrix(getSpace(), getOperator(common::Energy), getIndexConverter());
    }
    if (s_squared.has_value()) {
        s_squared->matrix_ =
            Matrix(getSpace(), getOperator(common::S_total_squared), getIndexConverter());
    }
    for (auto& [symbol_name, derivative] : derivative_of_energy_wrt_exchange_parameters) {
        // TODO: fix it!
        derivative.matrix_ = Matrix(
            getSpace(),
            getOperatorDerivative(common::Energy, model::symbols::J, symbol_name),
            getIndexConverter());
    }

    matrix_history_.matrices_was_built = true;
}

void Runner::BuildSpectra() {
    size_t number_of_blocks = getSpace().getBlocks().size();

    if (!getOperator(common::Energy).empty()) {
        energy.spectrum_.blocks.clear();
        energy.spectrum_.blocks.resize(number_of_blocks);
    }
    if (s_squared.has_value()) {
        s_squared->spectrum_.blocks.clear();
        s_squared->spectrum_.blocks.resize(number_of_blocks);
    }
    for (auto& [_, derivative] : derivative_of_energy_wrt_exchange_parameters) {
        derivative.spectrum_.blocks.clear();
        derivative.spectrum_.blocks.resize(number_of_blocks);
    }

    if (matrix_history_.matrices_was_built) {
        BuildSpectraUsingMatrices(number_of_blocks);
    } else {
        BuildSpectraWithoutMatrices(number_of_blocks);
    }
}

void Runner::BuildSpectraUsingMatrices(size_t number_of_blocks) {
    for (size_t block = 0; block < number_of_blocks; ++block) {
        DenseMatrix unitary_transformation_matrix;
        energy.spectrum_.blocks[block] =
            Subspectrum::energy(energy.matrix_.blocks[block], unitary_transformation_matrix);

        if (s_squared.has_value()) {
            s_squared->spectrum_.blocks[block] = Subspectrum::non_energy(
                s_squared->matrix_.blocks[block],
                unitary_transformation_matrix);
        }

        for (auto& [_, derivative] : derivative_of_energy_wrt_exchange_parameters) {
            derivative.spectrum_.blocks[block] = Subspectrum::non_energy(
                derivative.matrix_.blocks[block],
                unitary_transformation_matrix);
        }
    }
}

void Runner::BuildSpectraWithoutMatrices(size_t number_of_blocks) {
    for (size_t block = 0; block < number_of_blocks; ++block) {
        DenseMatrix unitary_transformation_matrix;
        {
            auto hamiltonian_submatrix = Submatrix(
                getSpace().getBlocks()[block],
                getOperator(common::Energy),
                getIndexConverter());
            energy.spectrum_.blocks[block] =
                Subspectrum::energy(hamiltonian_submatrix, unitary_transformation_matrix);
        }

        if (s_squared.has_value()) {
            auto non_hamiltonian_submatrix = Submatrix(
                getSpace().getBlocks()[block],
                getOperator(common::S_total_squared),
                getIndexConverter());
            s_squared->spectrum_.blocks[block] =
                Subspectrum::non_energy(non_hamiltonian_submatrix, unitary_transformation_matrix);
        }

        for (auto& [symbol_name, derivative] : derivative_of_energy_wrt_exchange_parameters) {
            // TODO: fix it
            auto derivative_submatrix = Submatrix(
                getSpace().getBlocks()[block],
                getOperatorDerivative(common::Energy, model::symbols::J, symbol_name),
                getIndexConverter());
            derivative.spectrum_.blocks[block] =
                Subspectrum::non_energy(derivative_submatrix, unitary_transformation_matrix);
        }
    }
}

const Matrix& Runner::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return energy.matrix_;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return s_squared->matrix_;
    }
}

const Spectrum& Runner::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return energy.spectrum_;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return s_squared->spectrum_;
    }
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
    DenseVector energy_vector;
    DenseVector degeneracy_vector;

    for (const auto& subspectrum : energy.spectrum_.blocks) {
        energy_vector.concatenate_with(subspectrum.raw_data);
        degeneracy_vector.add_identical_values(
            subspectrum.raw_data.size(),
            subspectrum.properties.degeneracy);
    }
    energy_vector.subtract_minimum();

    if (getSymbols().isAllGFactorsEqual()) {
        // and there is no field
        double g_factor = getSymbols().getGFactorParameters()->operator()(0);
        DenseVector s_squared_vector;

        // TODO: check if s_squared has been initialized
        for (const auto& subspectrum : s_squared->spectrum_.blocks) {
            s_squared_vector.concatenate_with(subspectrum.raw_data);
        }

        mu_squared_worker =
            std::make_unique<magnetic_susceptibility::UniqueGOnlySSquaredMuSquaredWorker>(
                std::move(energy_vector),
                std::move(degeneracy_vector),
                std::move(s_squared_vector),
                g_factor);
    } else {
        throw std::invalid_argument("Different g factors are not supported now.");
    }

    if (experimental_values_worker_.has_value()) {
        mu_squared_worker.value()->initializeExperimentalValues(
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

    if (mu_squared_worker.has_value()) {
        mu_squared_worker.value()->initializeExperimentalValues(
            experimental_values_worker_.value());
    }
}

std::map<model::symbols::SymbolName, double> Runner::calculateTotalDerivatives() {
    // TODO: ony s_squared-based calculation supported

    std::map<model::symbols::SymbolName, double> answer;

    for (const auto& changeable_symbol : getSymbols().getChangeableNames(model::symbols::J)) {
        DenseVector derivative_vector;
        for (const auto& subspectrum :
             getSpectrumDerivative(common::Energy, model::symbols::J, changeable_symbol).blocks) {
            derivative_vector.concatenate_with(subspectrum.raw_data);
        }
        double value = mu_squared_worker.value()->calculateTotalDerivative(
            model::symbols::J,
            std::move(derivative_vector));
        answer[changeable_symbol] = value;
        //        std::cout << "dR^2/d" << changeable_symbol << " = "
        //                  << value << std::endl;
    }

    if (!getSymbols().getChangeableNames(model::symbols::g_factor).empty()) {
        model::symbols::SymbolName g_name =
            getSymbols().getChangeableNames(model::symbols::g_factor)[0];
        double value =
            mu_squared_worker.value()->calculateTotalDerivative(model::symbols::g_factor);
        answer[g_name] = value;
        //        std::cout << "dR^2/d" << g_name << " = " << value << std::endl;
    }
    return answer;
}

void Runner::minimizeResidualError() {
    std::vector<model::symbols::SymbolName> changeable_names = getSymbols().getChangeableNames();
    std::vector<double> changeable_values;
    changeable_values.reserve(changeable_names.size());
    for (const model::symbols::SymbolName& name : changeable_names) {
        changeable_values.push_back(getSymbols().getValueOfName(name));
    }

    // Function-adapter for STLBFGS library.
    // It should calculate residual error and derivatives at a point of changeable_values.
    std::function<void(const std::vector<double>&, double&, std::vector<double>&)> func_grad_eval =
        [this, capture0 = std::cref(changeable_names)](
            const std::vector<double>& changeable_values,
            double& residual_error,
            std::vector<double>& gradient) {
            stepOfRegression(capture0, changeable_values, residual_error, gradient);
        };

    STLBFGS::Optimizer opt {func_grad_eval};
    opt.verbose = false;
    // Run calculation from initial guess. STLBFGS updates changeable_values every iteration.
    opt.run(changeable_values);

    for (size_t i = 0; i < changeable_names.size(); ++i) {
        //        std::cout << changeable_names[i] << ": " << changeable_values[i] << std::endl;
    }
}

void Runner::stepOfRegression(
    const std::vector<model::symbols::SymbolName>& changeable_names,
    const std::vector<double>& changeable_values,
    double& residual_error,
    std::vector<double>& gradient) {
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
    residual_error = mu_squared_worker.value()->calculateResidualError();

    // Calculate derivatives...
    std::map<model::symbols::SymbolName, double> map_gradient = calculateTotalDerivatives();
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        // ...and write it to external variable:
        gradient[i] = map_gradient[changeable_names[i]];
    }
}

const magnetic_susceptibility::MuSquaredWorker& Runner::getMuSquaredWorker() const {
    return *mu_squared_worker.value();
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
}  // namespace runner