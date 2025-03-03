#ifndef SPINNER_WIGNERECKARTHELPER_H
#define SPINNER_WIGNERECKARTHELPER_H

#include "src/common/index_converter/s_squared/Level.h"
#include "src/common/index_converter/s_squared/OrderOfSummation.h"
#include "src/spin_algebra/ClebshGordanCalculator.h"

namespace model::operators::ito {
class WignerEckartHelper {
  public:

	WignerEckartHelper() = delete;

	explicit WignerEckartHelper(
		std::shared_ptr<const index_converter::s_squared::OrderOfSummation> order_of_summation
	);

	double clebsh_gordan_coefficient(
		double l1,
		double l2,
		double l3,
		double m1,
		double m2) const;

	double local_product(
		const std::vector<spin_algebra::Multiplicity>& mults, 
		const std::vector<uint8_t>& ranks) const;

	double total_9j_coefficient(
		const index_converter::s_squared::Level& left,
		const index_converter::s_squared::Level& right,
		const std::vector<uint8_t>& ranks,
		double local_prod) const;

	void construct_overlapping_levels(const index_converter::s_squared::Level& level, 
		const std::vector<uint8_t>& ranks,
		std::vector<index_converter::s_squared::Level>& answer) const;

  private:
	spin_algebra::ClebshGordanCalculator clebshGordanCalculator_;
	std::shared_ptr<const index_converter::s_squared::OrderOfSummation> order_of_summation_;
};
} // namespace model::operators::ito

#endif // SPINNER_T20ONECENTERTERM_H
