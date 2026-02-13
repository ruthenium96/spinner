#include "src/entities/data_structures/AbstractFactories.h"
#include "src/entities/data_structures/eigen/EigenFactories.h"

namespace spinner::linear_algebra {
std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>
AbstractDenseTransformAndDiagonalizeFactory::defaultFactory() {
    auto answer = std::make_shared<EigenDenseTransformAndDiagonalizeFactory>();
    return answer;
}
} // namespace spinner::linear_algebra