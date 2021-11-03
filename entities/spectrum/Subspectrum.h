#ifndef JULY_SUBSPECTRUM_H
#define JULY_SUBSPECTRUM_H

#include <entities/data_structures/DenseMatrix.h>
#include "entities/BlockProperties.h"

struct Subspectrum {
    BlockProperties properties;
    DenseVector raw_data;

    friend std::ostream &operator<<(std::ostream &os, const Subspectrum &subspectrum);
};


#endif //JULY_SUBSPECTRUM_H
