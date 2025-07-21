#include "AbstractEigendecompositor.h"

#include <cassert>
#include <functional>
#include <stdexcept>
#include "src/common/Logger.h"
#include "src/common/OneOrMany.h"
#include "src/common/Quantity.h"

namespace {
template <typename T, typename U>
OneOrMany<U> getEntity(
    size_t number_of_subspaces, 
    const std::vector<OneOrMany<std::reference_wrapper<const T>>>& data) {
    if (std::all_of(
        data.cbegin(), 
        data.cend(), 
        holdsOne<std::reference_wrapper<const T>>
        )       
    ) {
        std::vector<std::reference_wrapper<const T>> answer;
        for (const auto& el: data) {
            answer.push_back(getOneRef(el));
        }
        return U(std::move(answer));
    } else {
        // TODO: check that number_of_many is the same for all Many.
        auto it = std::find_if(
            data.cbegin(),
            data.cend(),
            holdsMany<std::reference_wrapper<const T>>
        );
        size_t number_of_many = getManyRef(*it).size();
        std::vector<std::vector<std::reference_wrapper<const T>>> vector_of_answers(number_of_many);
        for (const auto& el: data) {
            if (holdsOne(el)) {
                auto reference = getOneRef(el);
                for (int j = 0; j < number_of_many; ++j) {
                    vector_of_answers[j].push_back(reference);
                }
            } else {
                auto vector_of_references = getManyRef(el);
                for (int j = 0; j < number_of_many; ++j) {
                    vector_of_answers[j].push_back(vector_of_references[j]);
                }
            }
        }
        std::vector<U> answer;
        for (int i = 0; i < vector_of_answers.size(); ++i) {
            answer.emplace_back(std::move(vector_of_answers[i]));
        }
        return answer;
    }
}

template <typename T, typename U>
std::optional<OneOrMany<U>> getEntity(
    size_t number_of_subspaces, 
    std::function<std::optional<OneOrMany<std::reference_wrapper<const T>>>(size_t)> getter) {
    std::vector<std::optional<OneOrMany<std::reference_wrapper<const T>>>> mb_data;
    for (int i = 0; i < number_of_subspaces; ++i) {
        mb_data.push_back(getter(i));
    }

    if (std::none_of(mb_data.cbegin(), mb_data.cend(), [](const auto& a){return a.has_value();})) {
        return std::nullopt;
    }
    if (!std::all_of(mb_data.cbegin(), mb_data.cend(), [](const auto& a){return a.has_value();})) {
        throw std::logic_error("Some Subspectra were not defined, while some were defined!");
    }

    std::vector<OneOrMany<std::reference_wrapper<const T>>> data;
    for (const auto& mb_el : mb_data) {
        data.push_back(mb_el.value());
    }
    return getEntity<T, U>(number_of_subspaces, data);
}

} // namespace

namespace eigendecompositor {

void AbstractEigendecompositor::BuildSpectra(
    const std::map<common::QuantityEnum, std::shared_ptr<const model::operators::Operator>>&
        operators,
    const std::map<
        std::pair<common::QuantityEnum, model::symbols::SymbolName>,
        std::shared_ptr<const model::operators::Operator>>& derivatives_operators,
    const space::Space& space) {
    buildSpectraWasCalled = true;
    {
        auto operators_to_calculate = operators;
        auto derivatives_operators_to_calculate = derivatives_operators;
        number_of_subspaces_ = space.getBlocks().size();
        initialize(
            operators_to_calculate,
            derivatives_operators_to_calculate,
            number_of_subspaces_);
        assert(operators_to_calculate.size() == 0);
        assert(derivatives_operators_to_calculate.size() == 0);
    }

    // todo: MemoryManager will have the special function: it will estimate required
    //  memory for each block and pick num_threads based on the sum of max values
    //  thus, openblas_num_threads can be set as floor(max_threads / num_threads)

    common::Logger::debug_msg("Eigendecomposition of...");
#pragma omp parallel for shared(space, operators, derivatives_operators) default(shared)
    for (size_t i = 0; i < space.getBlocks().size(); ++i) {
        common::Logger::debug_msg("block {} has started", i);
        const auto& subspace = space.getBlocks().at(i);
        BuildSubspectra(i, subspace);
        common::Logger::debug_msg("block {} is finished", i);
    }
    common::Logger::separate(2, common::debug);

    finalize();
}

bool AbstractEigendecompositor::BuildSpectraWasCalled() const {
    return buildSpectraWasCalled;
}

std::optional<OneOrMany<SpectrumRef>> AbstractEigendecompositor::getSpectrum(common::QuantityEnum quantity_enum) const {
    return getEntity<Subspectrum, SpectrumRef>(
        number_of_subspaces_,
        [this, quantity_enum](size_t number_of_block){
            return getSubspectrum(quantity_enum, number_of_block);
        });
}

std::optional<OneOrMany<MatrixRef>> AbstractEigendecompositor::getMatrix(common::QuantityEnum quantity_enum) const {
    return getEntity<Submatrix, MatrixRef>(
        number_of_subspaces_, 
        [this, quantity_enum](size_t number_of_block){
            return getSubmatrix(quantity_enum, number_of_block);
        }
    );
}

std::optional<OneOrMany<SpectrumRef>>
AbstractEigendecompositor::getSpectrumDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name) const {
    return getEntity<Subspectrum, SpectrumRef>(
        number_of_subspaces_, 
        [this, quantity_enum, symbol_name](size_t number_of_block){
            return getSubspectrumDerivative(quantity_enum, symbol_name, number_of_block);
        }
    );
}

std::optional<OneOrMany<SpectrumRef>>
AbstractEigendecompositor::getSpectrumDerivativeProduct(common::QuantityEnum quantity_enum, common::QuantityEnum quantity_enum_derivative, const model::symbols::SymbolName& symbol_name) const {
    return getEntity<Subspectrum, SpectrumRef>(
        number_of_subspaces_, 
        [this, quantity_enum, quantity_enum_derivative, symbol_name](size_t number_of_block){
            return getSubspectrumDerivativeProduct(quantity_enum, quantity_enum_derivative, symbol_name, number_of_block);
        }
    );
}

std::optional<OneOrMany<MatrixRef>>
AbstractEigendecompositor::getMatrixDerivative(common::QuantityEnum quantity_enum, const model::symbols::SymbolName& symbol_name) const {
    return getEntity<Submatrix, MatrixRef>(
        number_of_subspaces_, 
        [this, quantity_enum, symbol_name](size_t number_of_block){
            return getSubmatrixDerivative(quantity_enum, symbol_name, number_of_block);
        }
    );
}

size_t AbstractEigendecompositor::getSubspectrumSize(common::QuantityEnum quantity_enum, size_t number_of_block) const {
    auto data = getSubspectrum(quantity_enum, number_of_block).value();
    if (holdsOne(data)) {
        return getOneRef(data).get().raw_data->size();
    } else {
        return getManyRef(data)[0].get().raw_data->size();
    }
}

OneOrMany<std::vector<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>>
AbstractEigendecompositor::getWeightsOfAllStates() const {
    std::vector<OneOrMany<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>>
        weights;
    for (size_t i = 0; i < number_of_subspaces_; ++i) {
        weights.push_back(getWeightsOfBlockStates(i));
    }
    return getEntity<std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>, 
        std::vector<std::reference_wrapper<const std::unique_ptr<quantum::linear_algebra::AbstractDenseVector>>>>(
        number_of_subspaces_, weights);
}

}  // namespace eigendecompositor