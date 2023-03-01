#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/data_structures/std/StdFactories.h"
#include "src/entities/data_structures/std/StdSparseSemiunitaryMatrix.h"

namespace quantum::linear_algebra {

std::shared_ptr<AbstractSparseSemiunitaryFactory>
AbstractSparseSemiunitaryFactory::defaultSparseFactory() {
    return std::make_shared<StdSparseSemiunitaryFactory>();
}
}  // namespace quantum::linear_algebra