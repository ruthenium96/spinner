#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/data_structures/arma/ArmaFactories.h"

namespace quantum::linear_algebra {

std::shared_ptr<AbstractSparseSemiunitaryFactory>
AbstractSparseSemiunitaryFactory::defaultSparseFactory() {
    return std::make_shared<ArmaSparseSemiunitaryFactory>();
}
}  // namespace quantum::linear_algebra