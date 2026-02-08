#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/data_structures/arma/ArmaFactories.h"

namespace spinner::linear_algebra {

std::shared_ptr<AbstractSparseTransformFactory>
AbstractSparseTransformFactory::defaultSparseFactory() {
    return std::make_shared<ArmaSparseTransformFactory>();
}
} // namespace spinner::linear_algebra