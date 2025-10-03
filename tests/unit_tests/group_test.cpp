#include "src/group/Group.h"
#include <algorithm>

#include "gtest/gtest.h"
#include "magic_enum.hpp"

static std::vector<group::Group::GroupType> constuct_group_types() {
    std::vector<group::Group::GroupType> answer;
    answer.push_back(group::Group::S2);
    for (unsigned int order = 3; order < 50; order++) {
        answer.push_back({group::Group::Dihedral, order});
    }
    return answer;
}

static auto group_types = constuct_group_types();

TEST(group_info_tests, size_of_group_in_form_of_generators_equals_group_size) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        EXPECT_EQ(group_info.group_in_form_of_generators.size(), group_info.group_size)
            << "The size of group '" << group_type.type_enum
            << "' does not equal the size of group_in_form_of_generators";
    }
}

TEST(group_info_tests, size_of_element_of_group_in_form_of_generators_equals_number_of_generators) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        for (uint32_t i = 0; i < group_info.group_in_form_of_generators.size(); ++i) {
            EXPECT_EQ(
                group_info.group_in_form_of_generators[i].size(),
                group_info.number_of_generators)
                << i << " element of the group '" << group_type.type_enum
                << "' in the form of generators has wrong size";
        }
    }
}

TEST(group_info_tests, size_of_dimension_of_representation_equals_number_of_representations) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        EXPECT_EQ(
            group_info.number_of_representations,
            group_info.dimension_of_representation.size())
            << "In the group '" << group_type.type_enum
            << "' size of dimension_of_representation does not equal to number_of_representations";
    }
}

TEST(group_info_tests, size_of_orders_of_generators_equals_number_of_generators) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        EXPECT_EQ(group_info.number_of_generators, group_info.orders_of_generators.size())
            << "In the group '" << group_type.type_enum
            << "' size of orders_of_generators does not equal to number_of_generators";
    }
}

TEST(group_info_tests, size_of_orders_of_elements_equals_group_size) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        EXPECT_EQ(group_info.group_size, group_info.orders_of_elements.size())
            << "In the group '" << group_type.type_enum
            << "' size of orders_of_elements does not equal to group_size";
    }
}

TEST(
    group_info_tests,
    size_of_number_of_projectors_of_representation_equals_number_of_representations) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        EXPECT_EQ(
            group_info.number_of_representations,
            group_info.number_of_projectors_of_representation.size())
            << "In the group '" << group_type.type_enum
            << "' size of number_of_projectors_of_representation does not equal to number_of_representations";
    }
}

TEST(group_info_tests, sum_of_dimensions_of_representations_equals_group_size) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        uint32_t sum_of_squares = 0;
        for (auto d : group_info.dimension_of_representation) {
            sum_of_squares += d * d;
        }
        EXPECT_EQ(group_info.group_size, sum_of_squares)
            << "Sum of representations' dimensions of group '" << group_type.type_enum
            << "' does not equal to group_size";
    }
}

TEST(group_info_tests, size_of_coefficients_of_projectors_equals_number_of_representations) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        EXPECT_EQ(
            group_info.number_of_representations,
            group_info.coefficients_of_projectors.size())
            << "Vector size of coefficients_of_projectors in the group '" << group_type.type_enum
            << "' does not equal to number_of_representations";
    }
}

TEST(group_info_tests, number_of_projectors_of_representation_equals_actual_number_of_projectors) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        for (size_t repr = 0; repr < group_info.number_of_representations; ++repr) {
            EXPECT_EQ(
                group_info.number_of_projectors_of_representation[repr],
                group_info.coefficients_of_projectors[repr].size())
                << "number_of_projectors_of_representation '" << repr << "' in the group '"
                << group_type.type_enum << "' does not equal actual number of projectors";
        }
    }
}

TEST(group_info_tests, size_of_projectors_equals_group_size) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        for (size_t repr = 0; repr < group_info.number_of_representations; ++repr) {
            for (size_t k = 0; k < group_info.number_of_projectors_of_representation[repr]; ++k) {
                EXPECT_EQ(
                    group_info.group_size,
                    group_info.coefficients_of_projectors[repr][k].size())
                    << "In the group '" << group_type.type_enum << "', representation '" << repr
                    << "', projector '" << k << "', vector size does not equal to group_size";
            }
        }
    }
}

TEST(group_info_tests, descending_order_of_generators_order) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        size_t last_order_of_generator = UINT64_MAX;
        for (uint8_t generator = 0; generator < group_info.number_of_generators; ++generator) {
            size_t order_of_generator = 0;
            for (auto& v : group_info.group_in_form_of_generators) {
                order_of_generator = std::max(order_of_generator, v[generator]);
                EXPECT_GE(last_order_of_generator, order_of_generator)
                    << "In the group '" << group_type.type_enum << "' generator '" << generator
                    << "' has element order bigger than '" << (generator - 1)
                    << "', it can cause problems";
            }
        }
    }
}

TEST(group_info_tests, projectors_orthogonality) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        std::vector<std::vector<double>> projectors;
        for (uint8_t repr = 0; repr < group_info.number_of_representations; ++repr) {
            for (uint8_t k = 0; k < group_info.number_of_projectors_of_representation[repr]; ++k) {
                projectors.push_back(group_info.coefficients_of_projectors[repr][k]);
            }
        }

        for (uint32_t i = 0; i < group_info.group_size; ++i) {
            for (uint32_t j = i + 1; j < group_info.group_size; ++j) {
                double accumulator = 0;
                for (size_t g = 0; g < group_info.group_size; ++g) {
                    accumulator += projectors[i][g] * projectors[j][g];
                }
                EXPECT_NEAR(accumulator, 0.0, 1e-12)
                << "In the group '" << group_type.type_enum << "' projectors '"
                << i << "' and '" << j << "' are not orthogonal";
            }
        }
    }
}

TEST(group_info_tests, 00_corresponds_to_full_symmetric_representation) {
    for (auto& group_type : group_types) {
        const group::Group::AlgebraicProperties& group_info =
            group::Group::return_group_info_by_group_type(group_type);
        EXPECT_EQ(group_info.coefficients_of_projectors[0].size(), 1)
            << "In the group '" << group_type.type_enum << "' first representation has "
            << group_info.coefficients_of_projectors[0].size()
            << " projectors, but must have only one, because it should be full-symmetric representation";
        double value = group_info.coefficients_of_projectors[0][0][0];
        for (auto d : group_info.coefficients_of_projectors[0][0]) {
            EXPECT_EQ(d, value)
                << "In the group '" << group_type.type_enum
                << "' projector of the first representation has different coefficients, "
                << "this means it is not full-symmetric projector";
        }
    }
}

TEST(group_info_tests, size_of_multiplication_table) {
    const group::Group::AlgebraicProperties& group_info =
        group::Group::return_group_info_by_group_type(group::Group::S2);
    EXPECT_EQ(
        group_info.cayley_table.size(),
        group_info.number_of_representations * group_info.number_of_representations)
        << "In the group '" << group::Group::S2 << "' the size of multiplication table "
        << group_info.cayley_table.size()
        << " does not equal to square of number of representations "
        << group_info.number_of_representations * group_info.number_of_representations;
}

TEST(group_info_tests, full_symmetric_representation_as_identity_of_representation_multiplication) {
    const group::Group::AlgebraicProperties& group_info =
        group::Group::return_group_info_by_group_type(group::Group::S2);
    for (size_t i = 0; i < group_info.number_of_representations; ++i) {
        auto representation = std::set<uint8_t>();
        representation.insert(i);
        EXPECT_EQ(group_info.cayley_table.at({0, i}), representation)
            << "In the multiplication table of group '" << group::Group::S2
            << "' product of full-symmetric representation and " << i << " does not equal to "
            << i;
        EXPECT_EQ(group_info.cayley_table.at({i, 0}), representation)
            << "In the multiplication table of group '" << group::Group::S2 << "' product of" << i
            << " and full-symmetric representation "
            << " does not equal to " << i;
    }
}

TEST(
    group_info_tests,
    sizes_of_direct_product_representations_equal_to_product_of_dimension_of_representation) {
    const group::Group::AlgebraicProperties group_info =
        group::Group::return_group_info_by_group_type(group::Group::S2);
    for (size_t i = 0; i < group_info.number_of_representations; ++i) {
        for (size_t j = 0; j < group_info.number_of_representations; ++j) {
            auto product_of_sizes_of_representation = group_info.dimension_of_representation[i]
                * group_info.dimension_of_representation[j];
            int size_of_direct_product_representations = 0;
            for (auto k : group_info.cayley_table.at({i, j})) {
                size_of_direct_product_representations +=
                    group_info.dimension_of_representation[k];
            }
            EXPECT_EQ(
                product_of_sizes_of_representation,
                size_of_direct_product_representations)
                << "In the group '" << group::Group::S2
                << "' dimension of direct product of representations " << i << " * " << j
                << " = " << size_of_direct_product_representations
                << "does not equal to product of dimensions of these representations "
                << product_of_sizes_of_representation;
        }
    }
}

TEST(group_tests, throw_wrong_number_of_generators) {
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1}, {0, 0, 0}}),
        group::InitializationError);
    EXPECT_THROW(group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0}}), group::InitializationError);
}

TEST(group_tests, throw_wrong_size_of_generators) {
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 1, 4}}),
        group::InitializationError);
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2}}),
        group::InitializationError);
}

TEST(group_tests, throw_duplicate_of_number_in_generator) {
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 2}}),
        group::InitializationError);
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0}, {2, 2, 2}}),
        group::InitializationError);
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 0, 0}, {0, 2, 1}}),
        group::InitializationError);
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{0, 0, 0}, {0, 2, 1}}),
        group::InitializationError);
}

TEST(group_tests, throw_contains_wrong_number) {
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0}, {0, 2, 3}}),
        group::InitializationError);
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 3}, {0, 2, 1}}),
        group::InitializationError);
}

TEST(group_tests, throw_wrong_order_of_generators) {
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 3, 0}, {0, 2, 1, 3}}),
        group::InitializationError);
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0, 3}, {1, 2, 0, 3}}),
        group::InitializationError);
}

TEST(group_tests, throw_wrong_order_of_elements) {
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0, 3}, {0, 1, 3, 2}}),
        group::InitializationError);
    EXPECT_THROW(
        group::Group group({group::Group::Dihedral, 3}, {{1, 2, 0, 3}, {3, 1, 2, 0}}),
        group::InitializationError);
}

TEST(group_tests, construct_orbits_of_mults) {
    std::vector<std::pair<group::Group, std::vector<std::set<size_t>>>> cases = {
        {group::Group(group::Group::S2, {{1, 0, 3, 2, 4}}), {{0, 1}, {2, 3}, {4}}},
        {group::Group(group::Group::S2, {{4, 3, 2, 1, 0}}), {{0, 4}, {1, 3}, {2}}},
        {group::Group(group::Group::S2, {{4, 1, 2, 3, 0}}), {{0, 4}, {1}, {2}, {3}}},
        {group::Group({group::Group::Dihedral, 3}, {{1, 2, 0, 3}, {0, 2, 1, 3}}), {{0, 1, 2}, {3}}},
    };
    for (const auto& [group, right_answer] : cases) {
        auto orbits_of_mults = group.construct_orbits_of_mults();
        std::sort(orbits_of_mults.begin(), orbits_of_mults.end());
        EXPECT_EQ(orbits_of_mults, right_answer);
    }
}

TEST(group_tests, noncommutation_of_P2s_in_triangle) {
    auto group_one = group::Group(group::Group::S2, {{0, 2, 1}});
    auto group_two = group::Group(group::Group::S2, {{1, 0, 2}});
    EXPECT_FALSE(group_one.do_groups_commute(group_two));
    EXPECT_FALSE(group_two.do_groups_commute(group_one));
}

TEST(group_tests, commutation_of_P2s_in_rectangle) {
    auto group_one = group::Group(group::Group::S2, {{1, 0, 3, 2}});
    auto group_two = group::Group(group::Group::S2, {{2, 3, 0, 1}});
    EXPECT_TRUE(group_one.do_groups_commute(group_two));
    EXPECT_TRUE(group_two.do_groups_commute(group_one));
}

TEST(group_tests, commutation_of_P3s_in_toroid) {
    auto group_one = group::Group({group::Group::Dihedral, 3}, {
                                                        {
                                                            3, 4, 5,
                                                            6, 7, 8,
                                                            0, 1, 2,
                                                        },
                                                        {
                                                            3, 4, 5,
                                                            0, 1, 2,
                                                            6, 7, 8,
                                                        }
                                                    });
    auto group_two = group::Group({group::Group::Dihedral, 3}, {
                                                        {
                                                            1, 2, 0,
                                                            4, 5, 3,
                                                            7, 8, 6,
                                                        },
                                                        {
                                                            1, 0, 2,
                                                            4, 3, 5,
                                                            7, 6, 8,
                                                        }
                                                    });
    EXPECT_TRUE(group_one.do_groups_commute(group_two));
    EXPECT_TRUE(group_two.do_groups_commute(group_one));
}