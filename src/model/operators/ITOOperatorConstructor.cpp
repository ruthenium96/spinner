#include "ITOOperatorConstructor.h"
#include <cmath>
#include <memory>
#include <stdexcept>
#include "src/model/operators/terms/LocalSSquaredOneCenterTerm.h"
#include "src/model/operators/terms/ito/T00TwoCenterTerm.h"
#include "src/model/operators/terms/ito/T20OneCenterTerm.h"
#include "src/model/operators/terms/ito/T20TwoCenterTerm.h"
#include "src/model/operators/terms/ito/TotalSSquaredTerm.h"

namespace model::operators {
ITOOperatorConstructor::ITOOperatorConstructor(std::shared_ptr<index_converter::s_squared::IndexConverter> converter) :
	converter_(converter) {}

void ITOOperatorConstructor::emplaceIsotropicExchangeLike(
	std::shared_ptr<Operator> hamiltonian, 
	std::shared_ptr<const TwoDNumericalParameters<double>> parameters) const {
    hamiltonian->emplace_back(
        std::make_unique<const operators::ito::T00TwoCenterTerm>(
            converter_,
            parameters,
            2 * sqrt(3))
    );
}

void ITOOperatorConstructor::emplaceZeroFieldSplittingLike(
	std::shared_ptr<Operator> hamiltonian,
	std::shared_ptr<const OneDNumericalParameters<double>> parameters) const {
	throw std::invalid_argument("ITOOperatorConstructor::emplaceZeroFieldSplittingLike was not implemented yet");
}

std::shared_ptr<Operator> ITOOperatorConstructor::constructSSquared() const {
	auto s_squared_operator_ = std::make_shared<Operator>();
	s_squared_operator_->emplace_back(std::make_unique<const ito::TotalSSquaredTerm>(converter_));
    return s_squared_operator_;
}

std::shared_ptr<Operator> ITOOperatorConstructor::constructGSzSquaredLike(
	std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
	std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) const {
    // (\sum_i g_i S_i^z)^2 = \sum_i g_i^2 S_i^{z, 2} + \sum_i \sum_j g_i g_j S_i^z S_j^z
	auto g_sz_squared_operator_ = std::make_shared<Operator>();
    
    // \sum_i g_i^2 S_i^{z, 2} = \sum_i g_i^2 s_i(s_i+1) / 3 + ...
    g_sz_squared_operator_->emplace_back(
        std::make_unique<LocalSSquaredOneCenterTerm>(converter_, diagonal_parameters, 1.0 / 3.0)
    );
    // ... + \sum_i g_i^2 T_0^{(2)}(2|i) * sqrt(2/3)
    g_sz_squared_operator_->emplace_back(
        std::make_unique<const ito::T20OneCenterTerm>(converter_, diagonal_parameters, sqrt(2.0/3.0))
    );

    // \sum_i \sum_j g_i g_j S_i^z S_j^z = \sum_i \sum_j g_i g_j T_0^{(0)}(11|ij) * (-sqrt(1/3)) + ...
    g_sz_squared_operator_->emplace_back(
        std::make_unique<const ito::T00TwoCenterTerm>(converter_, nondiagonal_parameters, -2.0 * sqrt(1.0/3.0))
    );
    // ... + \sum_i \sum_j g_i g_j T_0^{(2)}(11|ij) * sqrt(2/3)
    g_sz_squared_operator_->emplace_back(
        std::make_unique<const ito::T20TwoCenterTerm>(converter_, nondiagonal_parameters, 2 * sqrt(2.0/3.0)));
    // this twos from summation in Submatrix: \sum_{a=1}^N \sum_{b=a+1}^N
    return g_sz_squared_operator_;
}

} // namespace model::operators