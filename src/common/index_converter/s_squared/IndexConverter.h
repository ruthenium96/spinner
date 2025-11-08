#ifndef SPINNER_S2INDEXCONVERTER_H
#define SPINNER_S2INDEXCONVERTER_H

#include "src/common/index_converter/AbstractIndexConverter.h"
#include "src/spin_algebra/Multiplicity.h"
#include "OrderOfSummation.h"
#include "Level.h"
#include <memory>
#include <optional>
#include <vector>

namespace index_converter::s_squared {

class IndexConverter : public AbstractIndexConverter {
public:
    IndexConverter(const std::vector<spin_algebra::Multiplicity>& mults,
    const std::shared_ptr<const OrderOfSummation>& order_of_summation);

    spin_algebra::Multiplicity convert_index_to_total_multiplicity(uint32_t index) const;

    std::pair<const Level&, uint8_t> convert_index_to_state(uint32_t index) const;
    std::optional<uint32_t> convert_state_to_index(const Level& level, uint8_t projection) const;

    uint8_t convert_index_to_tz_projection(uint32_t index) const override;

    std::shared_ptr<const OrderOfSummation> getOrderOfSummation() const;

private:
    std::shared_ptr<const OrderOfSummation> order_of_summation_;
    std::vector<std::vector<Level>> s_squared_levels_;
    std::vector<size_t> cumulative_sum_;
};
} // namespace index_converter::s_squared

#endif // SPINNER_S2INDEXCONVERTER_H