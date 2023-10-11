#include "HashmapFactories.h"

#include "EmhashSparseSemiunitaryMatrix.h"
#include "EmhashSparseSymmetricMatrix.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractSparseSemiunitaryMatrix>
EmhashSparseTransformFactory::createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) {
    auto answer = std::make_unique<EmhashSparseSemiunitaryMatrix>();
    answer->resize(cols, rows);
    return answer;
}

std::unique_ptr<AbstractSymmetricMatrix>
EmhashSparseTransformFactory::createSparseSymmetricMatrix(uint32_t size) {
    auto answer = std::make_unique<EmhashSparseSymmetricMatrix>();
    answer->resize(size);
    return answer;
}
}  // namespace quantum::linear_algebra