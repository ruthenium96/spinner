#include "src/linalg_structures/AbstractFactories.h"
#include "src/linalg_structures/arma/ArmaFactories.h"

namespace spinner::linalg_structures {
std::shared_ptr<AbstractDenseTransformAndDiagonalizeFactory>
AbstractDenseTransformAndDiagonalizeFactory::defaultFactory() {
    auto answer = std::make_shared<armadillo::ArmaDenseTransformAndDiagonalizeFactory>();
    return answer;
}
} // namespace spinner::linalg_structures