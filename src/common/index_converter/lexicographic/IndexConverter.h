#ifndef SPINNER_INDEXES_H
#define SPINNER_INDEXES_H

#include <cstdint>
#include <vector>

#include "src/spin_algebra/Multiplicity.h"
#include "src/common/index_converter/AbstractIndexConverter.h"

namespace index_converter::lexicographic {
class IndexConverter : public AbstractIndexConverter {
  public:
    explicit IndexConverter(std::vector<spin_algebra::Multiplicity> mults);
    std::vector<uint8_t> convert_lex_index_to_all_sz_projections(uint32_t lex) const;
    uint8_t convert_lex_index_to_one_sz_projection(uint32_t lex, uint32_t center) const;
    uint32_t convert_sz_projections_to_lex_index(const std::vector<uint8_t>& nzs) const;
    std::vector<IndexWithSign> convert_index_to_permutated_indexes(uint32_t index, 
                                                              const group::Group& group) const override;

    uint8_t convert_index_to_tz_projection(uint32_t index) const override;

    // NB: one should check projection value before ladding
    uint32_t ladder_projection(uint32_t lex, uint32_t center, int ladder) const;

  private:
    std::vector<uint32_t> cumulative_product_;
};
}  // namespace index_converter::lexicographic

#endif  // SPINNER_INDEXES_H
