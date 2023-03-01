#include "StdFactories.h"

#include "StdSparseSemiunitaryMatrix.h"
namespace quantum::linear_algebra {
std::unique_ptr<AbstractSparseSemiunitaryMatrix>
StdSparseSemiunitaryFactory::createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) {
    auto answer = std::make_unique<StdSparseSemiunitaryMatrix>();
    answer->resize(cols, rows);
    return answer;
}
}  // namespace quantum::linear_algebra