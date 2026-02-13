#include "AllSparseSemiunitaryMatrixFactories.h"

#include "src/linalg_structures/arma/ArmaFactories.h"
#include "src/linalg_structures/hashmap/HashmapFactories.h"

using namespace spinner::linalg_structures;

std::vector<std::shared_ptr<AbstractSparseTransformFactory>>
constructAllSparseSemiunitaryMatrixFactories() {
    return {
        std::make_shared<ArmaSparseTransformFactory>(),
        std::make_shared<EmhashSparseTransformFactory>()};
}
