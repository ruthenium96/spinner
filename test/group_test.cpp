#include <groups/Group.h>
#include "gtest/gtest.h"
#include "groups/Group_Info.h"

std::vector<group::GroupNames> group_names = {group::S2, group::S3};

TEST(group_info_tests, size_of_group_in_form_of_generators_equals_group_size) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        EXPECT_EQ(group_info.group_in_form_of_generators.size(), group_info.group_size)
        << "The size of group '" << group_name << "' does not equal the size of group_in_form_of_generators";
    }
}

TEST(group_info_tests, size_of_element_of_group_in_form_of_generators_equals_number_of_generators) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        for (uint32_t i = 0; i < group_info.group_in_form_of_generators.size(); ++i) {
            EXPECT_EQ(group_info.group_in_form_of_generators[i].size(), group_info.number_of_generators)
            << i << " element of the group '" << group_name << "' in the form of generators has wrong size";
        }
    }
}

TEST(group_info_tests, size_of_dimension_of_representation_equals_number_of_representations) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        EXPECT_EQ(group_info.number_of_representations, group_info.dimension_of_representation.size())
        << "In the group '" << group_name <<
        "' size of dimension_of_representation does not equal to number_of_representations";

    }
}

TEST(group_info_tests, size_of_number_of_projectors_of_representation_equals_number_of_representations) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        EXPECT_EQ(group_info.number_of_representations, group_info.number_of_projectors_of_representation.size())
        << "In the group '" << group_name <<
        "' size of number_of_projectors_of_representation does not equal to number_of_representations";

    }
}

TEST(group_info_tests, sum_of_dimensions_of_representations_equals_group_size) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        uint32_t sum_of_squares = 0;
        for (auto d : group_info.dimension_of_representation) {
            sum_of_squares += d*d;
        }
        EXPECT_EQ(group_info.group_size, sum_of_squares)
        << "Sum of representations' dimensions of group '" << group_name << "' does not equal to group_size";
    }
}

TEST(group_info_tests, size_of_coefficients_of_projectors_equals_number_of_representations) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        EXPECT_EQ(group_info.number_of_representations, group_info.coefficients_of_projectors.size())
        << "Vector size of coefficients_of_projectors in the group '" << group_name <<
        "' does not equal to number_of_representations";

    }
}

TEST(group_info_tests, number_of_projectors_of_representation_equals_actual_number_of_projectors) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        for (size_t repr = 0; repr < group_info.number_of_representations; ++repr) {
            EXPECT_EQ(group_info.number_of_projectors_of_representation[repr], group_info.coefficients_of_projectors[repr].size())
            << "number_of_projectors_of_representation '" << repr << "' in the group '" << group_name <<
            "' does not equal actual number of projectors";

        }
    }
}

TEST(group_info_tests, size_of_projectors_equals_group_size) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        for (size_t repr = 0; repr < group_info.number_of_representations; ++repr) {
            for (size_t k = 0; k < group_info.number_of_projectors_of_representation[repr]; ++k) {
                EXPECT_EQ(group_info.group_size, group_info.coefficients_of_projectors[repr][k].size())
                << "In the group '" << group_name << "', representation '" << repr <<
                "', projector '"<< k << "', vector size does not equal to group_size";
            }
        }
    }
}

TEST(group_info_tests, descending_order_of_generators_order) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        size_t last_order_of_generator = UINT64_MAX;
        for (uint8_t generator = 0; generator < group_info.number_of_generators; ++generator) {
            size_t order_of_generator = 0;
            for (auto& v : group_info.group_in_form_of_generators) {
                order_of_generator = std::max(order_of_generator, v[generator]);
                EXPECT_GE(last_order_of_generator, order_of_generator)
                << "In the group '" << group_name << "' generator '" << generator << "' has element order bigger than '"
                << (generator - 1) << "', it can cause problems";
            }
        }
    }
}

TEST(group_info_tests, projectors_orthogonality) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
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
                EXPECT_EQ(accumulator, 0.0) << "In the group '" << group_name << "' projectors '" << i << "' and '"
                << j << "' are not orthogonal";
            }
        }
    }
}

TEST(group_info_tests, 00_corresponds_to_full_symmetric_representation) {
    for (auto& group_name : group_names) {
        const group::GroupInfo& group_info = group::return_group_info_by_group_name(group_name);
        EXPECT_EQ(group_info.coefficients_of_projectors[0].size(), 1)
        << "In the group '" << group_name << "' first representation has "
        << group_info.coefficients_of_projectors[0].size()
        << " projectors, but must have only one, because it should be full-symmetric representation";
        double value = group_info.coefficients_of_projectors[0][0][0];
        for (auto d : group_info.coefficients_of_projectors[0][0]) {
            EXPECT_EQ(d, value) << "In the group '" << group_name
            << "' projector of the first representation has different coefficients, this means it is not full-symmetric projector";
        }
    }
}

TEST(group_tests, throw_wrong_number_of_generators) {
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0}, {0, 2, 1}, {0, 0, 0}}), std::length_error);
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0}}), std::length_error);
}

TEST(group_tests, throw_wrong_size_of_generators) {
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0}, {0, 2, 1, 4}}), std::length_error);
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0}, {0, 2}}), std::length_error);
}

TEST(group_tests, throw_duplicate_of_number_in_generator) {
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0}, {0, 2, 2}}), std::invalid_argument);
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0}, {2, 2, 2}}), std::invalid_argument);
    EXPECT_THROW(Group group(group::S3, {{1, 0, 0}, {0, 2, 1}}), std::invalid_argument);
    EXPECT_THROW(Group group(group::S3, {{0, 0, 0}, {0, 2, 1}}), std::invalid_argument);
}

TEST(group_tests, throw_contains_wrong_number) {
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0}, {0, 2, 3}}), std::invalid_argument);
    EXPECT_THROW(Group group(group::S3, {{1, 2, 3}, {0, 2, 1}}), std::invalid_argument);
}

TEST(group_tests, throw_wrong_order_of_generators) {
    EXPECT_THROW(Group group(group::S3, {{1, 2, 3, 0}, {0, 2, 1, 3}}), std::invalid_argument);
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0, 3}, {1, 2, 0, 3}}), std::invalid_argument);
}

TEST(group_tests, throw_wrong_order_of_elements) {
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0, 3}, {0, 1, 3, 2}}), std::invalid_argument);
    EXPECT_THROW(Group group(group::S3, {{1, 2, 0, 3}, {3, 1, 2, 0}}), std::invalid_argument);
}