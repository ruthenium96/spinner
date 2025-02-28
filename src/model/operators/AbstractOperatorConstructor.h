#ifndef SPINNER_ABSTRACTOPERATORCONSTRUCTOR_H
#define SPINNER_ABSTRACTOPERATORCONSTRUCTOR_H

#include <memory>
#include "Operator.h"

#include "src/model/NumericalParameters.h"

namespace model::operators {
class AbstractOperatorConstructor {
  public:
    virtual ~AbstractOperatorConstructor() = default;
    virtual void emplaceIsotropicExchangeLike(
        std::shared_ptr<Operator> hamiltonian, 
        std::shared_ptr<const TwoDNumericalParameters<double>> parameters) const = 0;

    virtual void emplaceZeroFieldSplittingLike(
        std::shared_ptr<Operator> hamiltonian,
        std::shared_ptr<const OneDNumericalParameters<double>> parameters) const = 0;

    virtual std::shared_ptr<Operator> constructSSquared() const = 0;
    virtual std::shared_ptr<Operator> constructMSquared() const = 0;
    
    virtual std::shared_ptr<Operator> constructGSzSquaredLike(
        std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
        std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) const = 0;

    virtual std::shared_ptr<Operator> constructGSquaredT00Like(
        std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
        std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) const = 0;
};
} // namespace model::operators

#endif // SPINNER_ABSTRACTOPERATORCONSTRUCTOR_H