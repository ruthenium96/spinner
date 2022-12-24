#ifndef SPINNER_STDFACTORIES_H
#define SPINNER_STDFACTORIES_H

#include "src/entities/data_structures/AbstractFactories.h"

namespace quantum::linear_algebra {
class StdSparseSemiunitaryFactory: public AbstractSparseSemiunitaryFactory {
    std::unique_ptr<AbstractSparseSemiunitaryMatrix>
    createSparseSemiunitaryMatrix(uint32_t cols, uint32_t rows) override;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_STDFACTORIES_H
