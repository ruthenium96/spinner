#include <memory>
#include <vector>
#include "gtest/gtest.h"

#include "src/common/index_converter/lexicographic/IndexConverter.h"
#include "src/common/index_converter/s_squared/IndexConverter.h"
#include "src/common/index_converter/s_squared/OrderOfSummation.h"
#include "src/group/Group.h"
#include "src/spin_algebra/Multiplicity.h"


TEST(converter_reversibility, lexicographic_2222_333_2345_44444) {
	std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults =
        {{2, 2, 2, 2}, {3, 3, 3}, {2, 3, 4, 5}, {4, 4, 4, 4, 4}};
	for (const auto& mults : vector_of_mults) {
		auto converter = std::make_shared<index_converter::lexicographic::IndexConverter>(mults);
		for (int index = 0; index < converter->get_total_space_size(); ++index) {
			auto state = converter->convert_lex_index_to_all_sz_projections(index);
			auto index_again = converter->convert_sz_projections_to_lex_index(state);
			EXPECT_EQ(index, index_again);
		}
	}
}

TEST(converter_reversibility, s_squared_2222_3333_4444) {
	std::vector<std::vector<spin_algebra::Multiplicity>> vector_of_mults =
        {{2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};
	std::vector<std::shared_ptr<const index_converter::s_squared::OrderOfSummation>> orders_of_summation;
	group::Group group_one(group::Group::S2, {{1, 0, 3, 2}});
	group::Group group_two(group::Group::S2, {{3, 2, 1, 0}});
	auto orbits_one = group_one.construct_orbits_of_mults();
	auto orbits_two = group_two.construct_orbits_of_mults();

	orders_of_summation.push_back(index_converter::s_squared::OrderOfSummation::constructFromOrbits(std::vector<std::vector<std::set<size_t>>>(), 4, 3));
	orders_of_summation.push_back(index_converter::s_squared::OrderOfSummation::constructFromOrbits(std::vector<std::vector<std::set<size_t>>>({orbits_one}), 4, 3));
	orders_of_summation.push_back(index_converter::s_squared::OrderOfSummation::constructFromOrbits(std::vector<std::vector<std::set<size_t>>>({orbits_two}), 4, 3));
	orders_of_summation.push_back(index_converter::s_squared::OrderOfSummation::constructFromOrbits(std::vector<std::vector<std::set<size_t>>>({orbits_one, orbits_two}), 4, 3));
	for (const auto& mults : vector_of_mults) {
		for (const auto& order_of_summation : orders_of_summation) {
		auto converter = std::make_shared<index_converter::s_squared::IndexConverter>(mults, order_of_summation);
			for (int index = 0; index < converter->get_total_space_size(); ++index) {
				auto state = converter->convert_index_to_state(index);
				auto index_again = converter->convert_state_to_index(state.first, state.second);
				EXPECT_TRUE(index_again.has_value());
				EXPECT_EQ(index, index_again.value());
			}
		}
	}
}

// TODO: the same test for s_squared converter