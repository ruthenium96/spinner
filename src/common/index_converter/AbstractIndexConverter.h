#ifndef SPINNER_ABSTRACT_INDEX_CONVERTER_H
#define SPINNER_ABSTRACT_INDEX_CONVERTER_H

#include <vector>
#include "src/spin_algebra/Multiplicity.h"
#include "src/group/Group.h"


namespace index_converter {

struct IndexWithSign {
    uint32_t index;
    int8_t sign;
};

class AbstractIndexConverter {
  public:
    explicit AbstractIndexConverter(std::vector<spin_algebra::Multiplicity> mults);
    const std::vector<spin_algebra::Multiplicity>& get_mults() const;
    const std::vector<double>& get_spins() const;
    uint32_t get_total_space_size() const;

    virtual uint8_t convert_index_to_tz_projection(uint32_t index) const = 0;
    virtual std::vector<IndexWithSign> convert_index_to_permutated_indexes(uint32_t index, 
                                                              const group::Group& group) const = 0;

    uint32_t get_max_ntz_proj() const;

  private:
	  std::vector<spin_algebra::Multiplicity> mults_;
    std::vector<double> spins_;
    uint32_t total_space_size_;
    uint32_t max_ntz_proj_;
};
} // namespace index_converter


#endif  // SPINNER_ABSTRACT_INDEX_CONVERTER_H