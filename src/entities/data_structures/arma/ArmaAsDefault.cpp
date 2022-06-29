#include "ArmaFactory.h"
#include "src/entities/data_structures/AbstractFactory.h"

namespace quantum::linear_algebra {
std::shared_ptr<AbstractFactory> AbstractFactory::defaultFactory() {
    auto answer = std::make_shared<ArmaFactory>();
    return answer;
}
}  // namespace quantum::linear_algebra