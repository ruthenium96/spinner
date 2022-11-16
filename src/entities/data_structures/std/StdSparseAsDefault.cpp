#include "StdSparseMatrix.h"
#include "src/entities/data_structures/AbstractSparseMatrix.h"

namespace quantum::linear_algebra {
std::unique_ptr<AbstractSparseMatrix> AbstractSparseMatrix::defaultSparseMatrix() {
    return std::make_unique<StdSparseMatrix>();
}
}  // namespace quantum::linear_algebra