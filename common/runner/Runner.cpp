#include "Runner.h"

#include <stlbfgs.h>

#include "components/matrix/MatrixBuilder.h"
#include "components/operator/ConstantOperator.h"
#include "components/operator/ScalarProduct.h"
#include "components/space/NonAbelianSimplifier.h"
#include "components/space/Symmetrizer.h"
#include "components/space/TzSorter.h"
#include "components/spectrum/SpectrumBuilder.h"

runner::Runner::Runner(const std::vector<int>& mults) :
    symbols_(mults.size()),
    converter_(mults),
    space_(converter_.get_total_space_size()) {
    operator_energy = Operator();
}

void runner::Runner::NonAbelianSimplify() {
    // TODO: model_is_finished
    if (space_history_.number_of_non_simplified_abelian_groups == 0) {
        return;
    }
    if (space_history_.number_of_non_simplified_abelian_groups != 1) {
        throw std::invalid_argument(
            "Non-Abelian simplification after using of two Non-Abelian Symmetrizers "
            "currently is not allowed.");
    }
    NonAbelianSimplifier nonAbelianSimplifier;
    space_ = nonAbelianSimplifier.apply(std::move(space_));
    space_history_.number_of_non_simplified_abelian_groups = 0;
    space_history_.isNonAbelianSimplified = true;
}

void runner::Runner::Symmetrize(group::Group new_group) {
    // TODO: model_is_finished
    // check if user trying to use the same Group for a second time:
    if (std::count(
            space_history_.applied_groups.begin(),
            space_history_.applied_groups.end(),
            new_group)) {
        return;
    }
    // TODO: symmetrizer does not work correct after non-Abelian simplifier. Fix it.
    if (space_history_.isNonAbelianSimplified && !new_group.properties.is_abelian) {
        throw std::invalid_argument(
            "Symmetrization after using of non-Abelian simplifier causes bugs.");
    }

    Symmetrizer symmetrizer(converter_, new_group);
    space_ = symmetrizer.apply(std::move(space_));

    if (!new_group.properties.is_abelian) {
        ++space_history_.number_of_non_simplified_abelian_groups;
    }
    space_history_.applied_groups.emplace_back(std::move(new_group));
}

void runner::Runner::Symmetrize(
    group::Group::GroupTypeEnum group_name,
    std::vector<group::Permutation> generators) {
    group::Group new_group(group_name, std::move(generators));
    Symmetrize(new_group);
}

void runner::Runner::TzSort() {
    // TODO: model_is_finished
    // It does not make any sense to use tz_sorter twice.
    if (space_history_.isTzSorted) {
        return;
    }
    TzSorter tz_sorter(converter_);
    space_ = tz_sorter.apply(std::move(space_));
    space_history_.isTzSorted = true;
}

const Space& runner::Runner::getSpace() const {
    return space_;
}

void runner::Runner::AssignSymbolToIsotropicExchange(
    const symbols::SymbolName& symbol_name,
    size_t center_a,
    size_t center_b) {
    if (model_is_finished) {
        throw std::invalid_argument("Adding parameters after finalization");
    }

    symbols_.assignSymbolToIsotropicExchange(symbol_name, center_a, center_b);
    if (!hamiltonian_history_.has_isotropic_exchange_interactions_initialized) {
        operator_energy.two_center_terms.emplace_back(
            std::make_unique<const ScalarProduct>(symbols_.getIsotropicExchangeParameters()));
        hamiltonian_history_.has_isotropic_exchange_interactions_initialized = true;
    }
}

void runner::Runner::BuildMatrices() {
    // TODO: refactor this:
    model_is_finished = true;
    for (const auto& applied_group : space_history_.applied_groups) {
        if (!symbols_.symmetry_consistence(applied_group)) {
            throw std::invalid_argument("Symbols do not match applied symmetries");
            // TODO: should we rename this exception?
        }
    }

    if (!space_history_.isNormalized) {
        for (auto& subspace : space_.blocks) {
            // TODO: maybe, we can implement normalize as Space method
            subspace.decomposition.normalize();
        }
        space_history_.isNormalized = true;
    }

    MatrixBuilder matrix_builder(converter_);
    if (!operator_energy.empty()) {
        matrix_energy.blocks.clear();
        matrix_energy = matrix_builder.apply(space_, operator_energy);
    }
    if (!operator_s_squared.empty()) {
        matrix_s_squared.blocks.clear();
        matrix_s_squared = matrix_builder.apply(space_, operator_s_squared);
    }
    for (const auto& [symbol, operator_derivative] :
         operator_derivative_of_energy_wrt_exchange_parameters) {
        matrix_derivative_of_energy_wrt_exchange_parameters[symbol].blocks.clear();
        matrix_derivative_of_energy_wrt_exchange_parameters[symbol] =
            matrix_builder.apply(space_, operator_derivative);
    }

    matrix_history_.matrices_was_built = true;
}

void runner::Runner::InitializeSSquared() {
    Operator s_squared_operator_;
    double sum_of_s_squared = 0;
    for (double spin : converter_.get_spins()) {
        sum_of_s_squared += spin * (spin + 1);
    }
    s_squared_operator_.zero_center_terms.emplace_back(
        std::make_unique<const ConstantOperator>(sum_of_s_squared));
    s_squared_operator_.two_center_terms.emplace_back(
        std::make_unique<const ScalarProduct>(converter_.get_spins().size()));

    operator_s_squared = std::move(s_squared_operator_);
}

void runner::Runner::InitializeIsotropicExchangeDerivatives() {
    // TODO: do nothing if this function is called twice
    for (const auto& symbol : symbols_.getChangeableNames(symbols::SymbolTypeEnum::J)) {
        Operator operator_derivative = Operator();
        operator_derivative.two_center_terms.emplace_back(std::make_unique<const ScalarProduct>(
            symbols_.constructIsotropicExchangeDerivativeParameters(symbol)));
        operator_derivative_of_energy_wrt_exchange_parameters[symbol] =
            std::move(operator_derivative);
    }
}

void runner::Runner::BuildSpectra() {
    model_is_finished = true;
    if (matrix_history_.matrices_was_built) {
        BuildSpectraUsingMatrices();
    } else {
        BuildSpectraWithoutMatrices();
    }
}

void runner::Runner::BuildSpectraUsingMatrices() {
    SpectrumBuilder spectrumBuilder;

    size_t number_of_blocks = space_.blocks.size();

    if (!operator_energy.empty()) {
        spectrum_energy.blocks.clear();
        spectrum_energy.blocks.resize(number_of_blocks);
    }
    if (!operator_s_squared.empty()) {
        spectrum_s_squared.blocks.clear();
        spectrum_s_squared.blocks.resize(number_of_blocks);
    }
    for (const auto& [symbol, _] : matrix_derivative_of_energy_wrt_exchange_parameters) {
        spectrum_derivative_of_energy_wrt_exchange_parameters[symbol] = Spectrum();
        spectrum_derivative_of_energy_wrt_exchange_parameters[symbol].blocks.resize(
            number_of_blocks);
    }

    for (size_t block = 0; block < number_of_blocks; ++block) {
        DenseMatrix unitary_transformation_matrix;
        spectrum_energy.blocks[block] = spectrumBuilder.apply_to_subentity_energy(
            matrix_energy.blocks[block],
            unitary_transformation_matrix);

        if (!operator_s_squared.empty()) {
            spectrum_s_squared.blocks[block] = spectrumBuilder.apply_to_subentity_non_energy(
                matrix_s_squared.blocks[block],
                unitary_transformation_matrix);
        }

        for (const auto& [symbol, matrix_derivative] :
             matrix_derivative_of_energy_wrt_exchange_parameters) {
            spectrum_derivative_of_energy_wrt_exchange_parameters[symbol].blocks[block] =
                spectrumBuilder.apply_to_subentity_non_energy(
                    matrix_derivative.blocks[block],
                    unitary_transformation_matrix);
        }
    }
}

void runner::Runner::BuildSpectraWithoutMatrices() {
    MatrixBuilder matrixBuilder(converter_);
    SpectrumBuilder spectrumBuilder;

    size_t number_of_blocks = space_.blocks.size();

    if (!operator_energy.empty()) {
        spectrum_energy.blocks.clear();
        spectrum_energy.blocks.resize(number_of_blocks);
    }
    if (!operator_s_squared.empty()) {
        spectrum_s_squared.blocks.clear();
        spectrum_s_squared.blocks.resize(number_of_blocks);
    }
    for (const auto& [symbol, _] : operator_derivative_of_energy_wrt_exchange_parameters) {
        spectrum_derivative_of_energy_wrt_exchange_parameters[symbol] = Spectrum();
        spectrum_derivative_of_energy_wrt_exchange_parameters[symbol].blocks.resize(
            number_of_blocks);
    }

    if (!space_history_.isNormalized) {
        for (auto& subspace : space_.blocks) {
            // TODO: maybe, we can implement normalize as Space method
            subspace.decomposition.normalize();
        }
        space_history_.isNormalized = true;
    }

    for (size_t block = 0; block < number_of_blocks; ++block) {
        DenseMatrix unitary_transformation_matrix;
        {
            Submatrix hamiltonian_submatrix =
                matrixBuilder.apply_to_subentity(space_.blocks[block], operator_energy);
            spectrum_energy.blocks[block] = spectrumBuilder.apply_to_subentity_energy(
                hamiltonian_submatrix,
                unitary_transformation_matrix);
        }

        if (!operator_s_squared.empty()) {
            Submatrix non_hamiltonian_submatrix =
                matrixBuilder.apply_to_subentity(space_.blocks[block], operator_s_squared);
            spectrum_s_squared.blocks[block] = spectrumBuilder.apply_to_subentity_non_energy(
                non_hamiltonian_submatrix,
                unitary_transformation_matrix);
        }

        for (const auto& [symbol, operator_derivative] :
             operator_derivative_of_energy_wrt_exchange_parameters) {
            Submatrix derivative_submatrix =
                matrixBuilder.apply_to_subentity(space_.blocks[block], operator_derivative);
            spectrum_derivative_of_energy_wrt_exchange_parameters[symbol].blocks[block] =
                spectrumBuilder.apply_to_subentity_non_energy(
                    derivative_submatrix,
                    unitary_transformation_matrix);
        }
    }
}

const Matrix& runner::Runner::getMatrix(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return matrix_energy;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return matrix_s_squared;
    }
}

const Spectrum& runner::Runner::getSpectrum(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return spectrum_energy;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return spectrum_s_squared;
    }
}

const Operator& runner::Runner::getOperator(common::QuantityEnum quantity_enum) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        return operator_energy;
    } else if (quantity_enum == common::QuantityEnum::S_total_squared) {
        return operator_s_squared;
    }
}

const lexicographic::IndexConverter& runner::Runner::getIndexConverter() const {
    return converter_;
}

void runner::Runner::AssignSymbolToGFactor(
    const symbols::SymbolName& symbol_name,
    size_t center_a) {
    // TODO: implement the same throws as in AssignSymbolToIsotropicExchange()
    symbols_.assignSymbolToGFactor(symbol_name, center_a);
}

const Operator& runner::Runner::getOperatorDerivative(
    common::QuantityEnum quantity_enum,
    symbols::SymbolTypeEnum symbol_type,
    const symbols::SymbolName& symbol) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        if (symbol_type == symbols::SymbolTypeEnum::J) {
            return operator_derivative_of_energy_wrt_exchange_parameters.at(symbol);
        }
    }
}

const Spectrum& runner::Runner::getSpectrumDerivative(
    common::QuantityEnum quantity_enum,
    symbols::SymbolTypeEnum symbol_type,
    const symbols::SymbolName& symbol) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        if (symbol_type == symbols::SymbolTypeEnum::J) {
            return spectrum_derivative_of_energy_wrt_exchange_parameters.at(symbol);
        }
    }
}

const Matrix& runner::Runner::getMatrixDerivative(
    common::QuantityEnum quantity_enum,
    symbols::SymbolTypeEnum symbol_type,
    const symbols::SymbolName& symbol) const {
    if (quantity_enum == common::QuantityEnum::Energy) {
        if (symbol_type == symbols::SymbolTypeEnum::J) {
            return matrix_derivative_of_energy_wrt_exchange_parameters.at(symbol);
        }
    }
}

void runner::Runner::BuildMuSquaredWorker() {
    DenseVector energy;
    DenseVector degeneracy;

    for (const auto& subspectrum : spectrum_energy.blocks) {
        energy.concatenate_with(subspectrum.raw_data);
        degeneracy.add_identical_values(
            subspectrum.raw_data.size(),
            subspectrum.properties.degeneracy);
    }
    energy.subtract_minimum();

    if (symbols_.isAllGFactorsEqual()) {
        // and there is no field
        double g_factor = symbols_.getGFactorParameters()->operator()(0);
        DenseVector s_squared;

        // TODO: check if s_squared has been initialized
        for (const auto& subspectrum : spectrum_s_squared.blocks) {
            s_squared.concatenate_with(subspectrum.raw_data);
        }

        mu_squared_worker =
            std::make_unique<magnetic_susceptibility::UniqueGOnlySSquaredMuSquaredWorker>(
                std::move(energy),
                std::move(degeneracy),
                std::move(s_squared),
                g_factor);
    } else {
        throw std::invalid_argument("Different g factors are not supported now.");
    }

    if (experimental_values_worker_.has_value()) {
        mu_squared_worker.value()->initializeExperimentalValues(
            experimental_values_worker_.value());
    }
}

void runner::Runner::initializeExperimentalValues(
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

std::map<symbols::SymbolName, double> runner::Runner::calculateTotalDerivatives() {
    // TODO: ony s_squared-based calculation supported

    std::map<symbols::SymbolName, double> answer;

    for (const auto& changeable_symbol : symbols_.getChangeableNames(symbols::J)) {
        DenseVector derivative;
        for (const auto& subspectrum :
             getSpectrumDerivative(common::Energy, symbols::J, changeable_symbol).blocks) {
            derivative.concatenate_with(subspectrum.raw_data);
        }
        double value =
            mu_squared_worker.value()->calculateTotalDerivative(symbols::J, std::move(derivative));
        answer[changeable_symbol] = value;
        //        std::cout << "dR^2/d" << changeable_symbol << " = "
        //                  << value << std::endl;
    }

    if (!symbols_.getChangeableNames(symbols::g_factor).empty()) {
        symbols::SymbolName g_name = symbols_.getChangeableNames(symbols::g_factor)[0];
        double value = mu_squared_worker.value()->calculateTotalDerivative(symbols::g_factor);
        answer[g_name] = value;
        //        std::cout << "dR^2/d" << g_name << " = " << value << std::endl;
    }
    return answer;
}

void runner::Runner::minimizeResidualError() {
    std::vector<symbols::SymbolName> changeable_names = symbols_.getChangeableNames();
    std::vector<double> changeable_values;
    changeable_values.reserve(changeable_names.size());
    for (const symbols::SymbolName& name : changeable_names) {
        changeable_values.push_back(symbols_.getValueOfName(name));
    }

    std::function<void(const std::vector<double>&, double&, std::vector<double>&)> func_grad_eval =
        [this, &changeable_names](
            const std::vector<double>& changeable_values,
            double& residual_error,
            std::vector<double>& gradient) {
            for (size_t i = 0; i < changeable_names.size(); ++i) {
                symbols_.setNewValueToChangeableSymbol(changeable_names[i], changeable_values[i]);
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

            std::map<symbols::SymbolName, double> map_gradient = calculateTotalDerivatives();
            for (size_t i = 0; i < changeable_names.size(); ++i) {
                gradient[i] = map_gradient[changeable_names[i]];
            }
        };

    STLBFGS::Optimizer opt {func_grad_eval};
    opt.verbose = false;
    opt.run(changeable_values);

    for (size_t i = 0; i < changeable_names.size(); ++i) {
        //        std::cout << changeable_names[i] << ": " << changeable_values[i] << std::endl;
    }
}

const std::unique_ptr<magnetic_susceptibility::MuSquaredWorker>&
runner::Runner::getPtrToMuSquaredWorker() const {
    return mu_squared_worker.value();
}

const symbols::Symbols& runner::Runner::getSymbols() const {
    return symbols_;
}
symbols::Symbols& runner::Runner::modifySymbols() {
    // TODO: model_is_finished
    return symbols_;
}
