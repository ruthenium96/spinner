#include "AllSparseSemiunitaryMatrixFactories.h"

#include "src/entities/data_structures/arma/ArmaFactories.h"
#include "src/entities/data_structures/hashmap/HashmapFactories.h"
std::vector<std::shared_ptr<quantum::linear_algebra::AbstractSparseTransformFactory>>
constructAllSparseSemiunitaryMatrixFactories() {
    return {
        std::make_shared<quantum::linear_algebra::ArmaSparseTransformFactory>(),
        std::make_shared<quantum::linear_algebra::EmhashSparseTransformFactory>()};
}
