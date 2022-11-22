#include "EigenDenseFactory.h"
#include "src/entities/data_structures/AbstractDenseFactory.h"

namespace quantum::linear_algebra {
std::shared_ptr<AbstractDenseFactory> AbstractDenseFactory::defaultFactory() {
    auto answer = std::make_shared<EigenDenseFactory>();
    return answer;
}
}  // namespace quantum::linear_algebra