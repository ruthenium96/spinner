#ifndef SPINNER_EMHASHSPARSESYMMETRICMATRIX_H
#define SPINNER_EMHASHSPARSESYMMETRICMATRIX_H

#include "hash_table7.hpp"
#include "src/entities/data_structures/AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {

class EmhashSparseSymmetricMatrix: public AbstractSymmetricMatrix {
  public:
    void add_to_position(double value, uint32_t i, uint32_t j) override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const noexcept override;
    void print(std::ostream& os) const override;
    ~EmhashSparseSymmetricMatrix() override = default;
    void resize(size_t size);

  private:
    size_t size_;
    emhash7::HashMap<uint32_t, emhash7::HashMap<uint32_t, double>> hashmap_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EMHASHSPARSESYMMETRICMATRIX_H