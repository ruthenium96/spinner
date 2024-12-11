#ifndef SPINNER_TOTALSSQUAREDTERM_H
#define SPINNER_TOTALSSQUAREDTERM_H

#include "src/model/operators/terms/Term.h"
#include "src/common/index_converter/s_squared/IndexConverter.h"


namespace model::operators::ito {
class TotalSSquaredTerm : public ZeroCenterTerm {
  public:
    TotalSSquaredTerm(std::shared_ptr<const index_converter::s_squared::IndexConverter> converter);
    std::unique_ptr<Term> clone() const override;
    void construct(
        quantum::linear_algebra::AbstractSymmetricMatrix& matrix,
        const std::set<unsigned int>& indexes_of_vectors) const override;

  private:
        std::shared_ptr<const index_converter::s_squared::IndexConverter> converter_;

};
} // namespace model::operators::ito

#endif // SPINNER_TOTALSSQUAREDTERM_H