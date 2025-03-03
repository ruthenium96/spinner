#ifndef SPINNER_ITOOPERATORCONSTRUCTOR_H
#define SPINNER_ITOOPERATORCONSTRUCTOR_H

#include "AbstractOperatorConstructor.h"

#include "src/common/index_converter/s_squared/IndexConverter.h"

namespace model::operators {
class ITOOperatorConstructor : public AbstractOperatorConstructor {
  public:
    ITOOperatorConstructor(std::shared_ptr<index_converter::s_squared::IndexConverter> converter);

    void emplaceIsotropicExchangeLike(
        std::shared_ptr<Operator> hamiltonian, 
        std::shared_ptr<const TwoDNumericalParameters<double>> parameters) const override;

    void emplaceZeroFieldSplittingLike(
        std::shared_ptr<Operator> hamiltonian,
        std::shared_ptr<const OneDNumericalParameters<double>> parameters) const override;

    std::shared_ptr<Operator> constructSSquared() const override;
    std::shared_ptr<Operator> constructMSquared() const override;

    std::shared_ptr<Operator> constructGSzSquaredLike(
        std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
        std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) const override;

    std::shared_ptr<Operator> constructGSquaredT00Like(
        std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
        std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) const override;
  
  private:
    std::shared_ptr<index_converter::s_squared::IndexConverter> converter_;

};
} // namespace model::operators

#endif // SPINNER_ITOOPERATORCONSTRUCTOR_H