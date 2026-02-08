#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/data_structures/hashmap/HashmapFactories.h"

namespace spinner::quantum::linear_algebra {

std::shared_ptr<AbstractSparseTransformFactory>
AbstractSparseTransformFactory::defaultSparseFactory() {
    return std::make_shared<EmhashSparseTransformFactory>();
}
} // namespace spinner::quantum::linear_algebra