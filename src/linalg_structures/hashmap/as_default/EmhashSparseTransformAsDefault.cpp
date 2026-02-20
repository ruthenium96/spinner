#include "src/linalg_structures/AbstractFactories.h"
#include "src/linalg_structures/hashmap/HashmapFactories.h"

namespace spinner::linalg_structures {

std::shared_ptr<AbstractSparseTransformFactory>
AbstractSparseTransformFactory::defaultSparseFactory() {
    return std::make_shared<hashmap::EmhashSparseTransformFactory>();
}
} // namespace spinner::linalg_structures