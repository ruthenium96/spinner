#include "AllSparseSemiunitaryMatrixFactories.h"

#include "src/entities/data_structures/arma/ArmaFactories.h"
#include "src/entities/data_structures/hashmap/HashmapFactories.h"

using namespace spinner::linear_algebra;

std::vector<std::shared_ptr<AbstractSparseTransformFactory>>
constructAllSparseSemiunitaryMatrixFactories() {
    return {
        std::make_shared<ArmaSparseTransformFactory>(),
        std::make_shared<EmhashSparseTransformFactory>()};
}
