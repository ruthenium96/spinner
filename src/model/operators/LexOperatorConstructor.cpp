#include "LexOperatorConstructor.h"
#include <memory>
#include "src/model/operators/Operator.h"
#include "src/model/operators/terms/ConstantTerm.h"
#include "src/model/operators/terms/LocalSSquaredOneCenterTerm.h"
#include "src/model/operators/terms/lexicographic/ScalarProductTerm.h"
#include "src/model/operators/terms/lexicographic/SzSzOneCenterTerm.h"
#include "src/model/operators/terms/lexicographic/SzSzTwoCenterTerm.h"


namespace model::operators {
LexOperatorConstructor::LexOperatorConstructor(std::shared_ptr<index_converter::lexicographic::IndexConverter> converter) :
	converter_(converter) {}

void LexOperatorConstructor::emplaceIsotropicExchangeLike(
	std::shared_ptr<Operator> hamiltonian, 
	std::shared_ptr<const TwoDNumericalParameters<double>> parameters) const {
	hamiltonian->emplace_back(
	std::make_unique<const operators::lexicographic::ScalarProductTerm>(
		converter_,
		parameters));
}

void LexOperatorConstructor::emplaceZeroFieldSplittingLike(
	std::shared_ptr<Operator> hamiltonian,
	std::shared_ptr<const OneDNumericalParameters<double>> parameters) const {
	hamiltonian->emplace_back(
        std::make_unique<const operators::LocalSSquaredOneCenterTerm>(
            converter_,
            parameters,
            -1.0 / 3.0));
    hamiltonian->emplace_back(
        std::make_unique<const operators::lexicographic::SzSzOneCenterTerm>(converter_, parameters));
}

std::shared_ptr<Operator> LexOperatorConstructor::constructSSquared() const {
	auto s_squared_operator_ = std::make_shared<Operator>();
    auto sum_of_s_squared = std::make_shared<double>(0);
    for (double spin : converter_->get_spins()) {
        *sum_of_s_squared += spin * (spin + 1);
    }
	// (\sum_i S_i)^2 = \sum_i S_i^2 + ...
	s_squared_operator_->emplace_back(std::make_unique<const ConstantTerm>(sum_of_s_squared));
	// ... + 2 * \sum_i\sum_j (S_i, S_j). This ScalarProductTerm constructor does required job:
	s_squared_operator_->emplace_back(std::make_unique<const lexicographic::ScalarProductTerm>(converter_));
    return s_squared_operator_;
}

std::shared_ptr<Operator> LexOperatorConstructor::constructGSzSquaredLike(
	std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
	std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) const {
	auto g_sz_squared_operator_ = std::make_shared<Operator>();

    g_sz_squared_operator_->emplace_back(
        std::make_unique<const lexicographic::SzSzOneCenterTerm>(converter_, diagonal_parameters));
    g_sz_squared_operator_->emplace_back(
        std::make_unique<const lexicographic::SzSzTwoCenterTerm>(converter_, nondiagonal_parameters, 2));
    // this two from summation in Submatrix: \sum_{a=1}^N \sum_{b=a+1}^N
    return g_sz_squared_operator_;
}

} // namespace model::operators