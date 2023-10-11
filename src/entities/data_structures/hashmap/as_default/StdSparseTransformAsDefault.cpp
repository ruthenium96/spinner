#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/data_structures/hashmap/StdFactories.h"

namespace quantum::linear_algebra {

std::shared_ptr<AbstractSparseTransformFactory>
AbstractSparseTransformFactory::defaultSparseFactory() {
    return std::make_shared<StdSparseTransformFactory>();
}
}  // namespace quantum::linear_algebra