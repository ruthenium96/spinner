#ifndef JULY_SUBSPECTRUM_H
#define JULY_SUBSPECTRUM_H

#include "entities/BlockProperties.h"
#include "entities/data_structures/DenseVector.h"

struct Subspectrum {
    BlockProperties properties;
    DenseVector energy_raw_data;
    std::vector<DenseVector> non_energy_raw_data;

    friend std::ostream &operator<<(std::ostream &os, const Subspectrum &subspectrum);
};


#endif //JULY_SUBSPECTRUM_H
