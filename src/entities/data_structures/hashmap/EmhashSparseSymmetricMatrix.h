#ifndef SPINNER_EMHASHSPARSESYMMETRICMATRIX_H
#define SPINNER_EMHASHSPARSESYMMETRICMATRIX_H

#include "hash_table8.hpp"
#include "src/entities/data_structures/AbstractSymmetricMatrix.h"

namespace quantum::linear_algebra {

class EmhashSparseSymmetricMatrix: public AbstractSymmetricMatrix {
  public:
    using Map = emhash8::HashMap<uint32_t, double>;

    void add_to_position(double value, uint32_t i, uint32_t j) override;
    uint32_t size() const override;
    double at(uint32_t i, uint32_t j) const noexcept override;
    void print(std::ostream& os) const override;
    ~EmhashSparseSymmetricMatrix() override = default;
    void resize(size_t size);

    const emhash8::HashMap<uint32_t, Map>& getSparseSymmetricMatrix() const;

  private:
    size_t size_;
    emhash8::HashMap<uint32_t, Map> hashmap_;
};
}  // namespace quantum::linear_algebra
#endif  //SPINNER_EMHASHSPARSESYMMETRICMATRIX_H
