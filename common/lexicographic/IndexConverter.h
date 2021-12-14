#ifndef JULY_INDEXES_H
#define JULY_INDEXES_H

#include <cstdint>
#include <functional>
#include <numeric>
#include <vector>

namespace lexicographic {
class IndexConverter {
  public:
    explicit IndexConverter(std::vector<int> mults);
    [[nodiscard]] uint8_t convert_lex_index_to_tz_projection(uint32_t lex) const;
    [[nodiscard]] std::vector<uint8_t> convert_lex_index_to_all_sz_projections(uint32_t lex) const;
    [[nodiscard]] uint8_t
    convert_lex_index_to_one_sz_projection(uint32_t lex, uint32_t center) const;
    [[nodiscard]] uint32_t
    convert_sz_projections_to_lex_index(const std::vector<uint8_t>& nzs) const;

    // NB: one should check projection value before ladding
    [[nodiscard]] uint32_t ladder_projection(uint32_t lex, uint32_t center, int ladder) const;

    [[nodiscard]] const std::vector<int>& get_mults() const;
    [[nodiscard]] const std::vector<double>& get_spins() const;
    [[nodiscard]] uint32_t get_total_space_size() const;

  private:
    std::vector<uint32_t> cumulative_product;
    std::vector<int> mults_;
    std::vector<double> spins_;
    uint32_t total_space_size;
};
}  // namespace lexicographic

#endif  // JULY_INDEXES_H
