#ifndef SPINNER_LEXOPERATORCONSTRUCTOR_H
#define SPINNER_LEXOPERATORCONSTRUCTOR_H

#include <memory>
#include "AbstractOperatorConstructor.h"
#include "src/common/index_converter/lexicographic/IndexConverter.h"

namespace model::operators {
class LexOperatorConstructor : public AbstractOperatorConstructor {
  public:
    LexOperatorConstructor(std::shared_ptr<index_converter::lexicographic::IndexConverter> converter);

    void emplaceIsotropicExchangeLike(
        std::shared_ptr<Operator> hamiltonian, 
        std::shared_ptr<const TwoDNumericalParameters<double>> parameters) const override;

    void emplaceZeroFieldSplittingLike(
        std::shared_ptr<Operator> hamiltonian,
        std::shared_ptr<const OneDNumericalParameters<double>> parameters) const override;

    std::shared_ptr<Operator> constructSSquared() const override;
    std::shared_ptr<Operator> constructGSzSquaredLike(
        std::shared_ptr<const OneDNumericalParameters<double>> diagonal_parameters,
        std::shared_ptr<const TwoDNumericalParameters<double>> nondiagonal_parameters) const override;

  private:
    std::shared_ptr<index_converter::lexicographic::IndexConverter> converter_;

};
} // namespace model::operators

#endif // SPINNER_LEXOPERATORCONSTRUCTOR_H