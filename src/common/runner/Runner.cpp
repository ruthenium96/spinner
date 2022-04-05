#include "Runner.h"

#include <stlbfgs.h>

#include <utility>

#include "src/model/operators/ConstantTerm.h"
#include "src/model/operators/ScalarProductTerm.h"
#include "src/space/optimization/NonAbelianSimplifier.h"
#include "src/space/optimization/PositiveProjectionsEliminator.h"
#include "src/space/optimization/Symmetrizer.h"
#include "src/space/optimization/TzSorter.h"

namespace runner {

Runner::Runner(model::Model model) :
    Runner(std::move(model), common::physical_optimization::OptimizationList()) {}

Runner::Runner(
    model::Model model,
    common::physical_optimization::OptimizationList optimizationList) :
    model_(std::move(model)),
    optimizationList_(std::move(optimizationList)),
    // TODO: Move it from here, because we do not need to allocate Space if Model
    //  and OptimizationList are not consistence.
    space_(model_.getIndexConverter().get_total_space_size()) {
    if (model_.is_s_squared_initialized()) {
        s_squared = common::Quantity();
    }
    if (model_.is_isotropic_exchange_derivatives_initialized()) {
        for (const auto& symbol :
             getSymbols().getChangeableNames(model::symbols::SymbolTypeEnum::J)) {
            derivative_of_energy_wrt_exchange_parameters[symbol] = common::Quantity();
        }
    }
    // TODO: check here model_ and optimizationList_ consistence!

    for (const auto& applied_group : optimizationList_.getGroupsToApply()) {
        if (!getSymbols().symmetry_consistence(applied_group)) {
            throw std::invalid_argument("Symbols do not match applied symmetries");
            // TODO: should we rename this exception?
        }
    }

    if (getSymbols().isIsotropicExchangeInitialized()) {
        model_.InitializeIsotropicExchange();
    }

    //    if (!symbols_.isGFactorInitialized()) {
    //        throw std::length_error("g factor parameters have not been initialized");
    //    }

    if (optimizationList_.isTzSorted()) {
        TzSort();
    }
    if (optimizationList_.isPositiveProjectionsEliminated()) {
        EliminatePositiveProjections();
    }
    for (const auto& group : optimizationList_.getGroupsToApply()) {
        Symmetrize(group);
    }
    for (auto& subspace : space_.getBlocks()) {
        // TODO: maybe, we can implement normalize as Space method
        subspace.decomposition.normalize();
    }
    space_history_.isNormalized = true;
}

void Runner::EliminatePositiveProjections() {
    uint32_t max_ntz_proj = getIndexConverter().get_max_ntz_proj();

    space::optimization::PositiveProjectionsEliminator positiveProjectionsEliminator(max_ntz_proj);
    space_ = positiveProjectionsEliminator.apply(std::move(space_));
}

void Runner::NonAbelianSimplify() {
    if (space_history_.number_of_non_simplified_abelian_groups == 0) {
        return;
    }
    if (space_history_.number_of_non_simplified_abelian_groups != 1) {
        throw std::invalid_argument(
            "Non-Abelian simplification after using of two Non-Abelian Symmetrizers "
            "currently is not allowed.");
    }
    space::optimization::NonAbelianSimplifier nonAbelianSimplifier;
    space_ = nonAbelianSimplifier.apply(std::move(space_));
    space_history_.number_of_non_simplified_abelian_groups = 0;
    space_history_.isNonAbelianSimplified = true;
}

void Runner::Symmetrize(group::Group new_group) {
    // TODO: symmetrizer does not work correct after non-Abelian simplifier. Fix it.
    if (space_history_.isNonAbelianSimplified && !new_group.properties.is_abelian) {
        throw std::invalid_argument(
            "Symmetrization after using of non-Abelian simplifier causes bugs.");
    }

    space::optimization::Symmetrizer symmetrizer(getIndexConverter(), new_group);
    space_ = symmetrizer.apply(std::move(space_));

    if (!new_group.properties.is_abelian) {
        ++space_history_.number_of_non_simplified_abelian_groups;
    }
}

void Runner::TzSort() {
    space::optimization::TzSorter tz_sorter(getIndexConverter());
    space_ = tz_sorter.apply(std::move(space_));
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
    return model_.getOperator(quantity_enum);
}

const lexicographic::IndexConverter& Runner::getIndexConverter() const {
    return model_.getIndexConverter();
}

const model::operators::Operator& Runner::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    model::symbols::SymbolTypeEnum symbol_type,
    const model::symbols::SymbolName& symbol) const {
    return model_.getOperatorDerivative(quantity_enum, symbol_type, symbol);
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

    using namespace std::placeholders;
    std::function<void(const std::vector<double>&, double&, std::vector<double>&)> func_grad_eval =
        std::bind(&Runner::stepOfRegression, this, std::cref(changeable_names), _1, _2, _3);

    STLBFGS::Optimizer opt {func_grad_eval};
    opt.verbose = false;
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
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        // TODO: mutable use of Model/Symbols. Refactor it.
        model_.getSymbols().setNewValueToChangeableSymbol(
            changeable_names[i],
            changeable_values[i]);
    }

    if (matrix_history_.matrices_was_built) {
        // TODO: 1) it does not work now
        //       2) here's the problem with future SSquaredTransformer
        //       or we can use SSquaredTransform() as flag and actually apply it at BuildMatrix?
        BuildMatrices();
    }
    BuildSpectra();

    BuildMuSquaredWorker();

    residual_error = mu_squared_worker.value()->calculateResidualError();

    std::map<model::symbols::SymbolName, double> map_gradient = calculateTotalDerivatives();
    for (size_t i = 0; i < changeable_names.size(); ++i) {
        gradient[i] = map_gradient[changeable_names[i]];
    }
}

const magnetic_susceptibility::MuSquaredWorker& Runner::getMuSquaredWorker() const {
    return *mu_squared_worker.value();
}

const model::symbols::Symbols& Runner::getSymbols() const {
    return model_.getSymbols();
}
}  // namespace runner