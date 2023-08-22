#ifndef SPINNER_INDEXES_H
#define SPINNER_INDEXES_H

#include <cstdint>
#include <functional>
#include <numeric>
#include <vector>

#include "src/spin_algebra/Multiplicity.h"

namespace lexicographic {
class IndexConverter {
  public:
    explicit IndexConverter(std::vector<spin_algebra::Multiplicity> mults);
    uint8_t convert_lex_index_to_tz_projection(uint32_t lex) const;
    std::vector<uint8_t> convert_lex_index_to_all_sz_projections(uint32_t lex) const;
    uint8_t convert_lex_index_to_one_sz_projection(uint32_t lex, uint32_t center) const;
    uint32_t convert_sz_projections_to_lex_index(const std::vector<uint8_t>& nzs) const;

    // NB: one should check projection value before ladding
    uint32_t ladder_projection(uint32_t lex, uint32_t center, int ladder) const;

    const std::vector<spin_algebra::Multiplicity>& get_mults() const;
    const std::vector<double>& get_spins() const;
    uint32_t get_total_space_size() const;
    uint32_t get_max_ntz_proj() const;

  private:
    std::vector<uint32_t> cumulative_product;
    std::vector<spin_algebra::Multiplicity> mults_;
    std::vector<double> spins_;
    uint32_t total_space_size;
};
}  // namespace lexicographic

#endif  // SPINNER_INDEXES_H
