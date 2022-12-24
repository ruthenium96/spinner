#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/data_structures/arma/ArmaFactories.h"

namespace quantum::linear_algebra {
std::shared_ptr<AbstractDenseVectorFactory> AbstractDenseVectorFactory::defaultFactory() {
    auto answer = std::make_shared<ArmaDenseVectorFactory>();
    return answer;
}
}  // namespace quantum::linear_algebra