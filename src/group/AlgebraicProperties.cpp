#include <cmath>
#include <cstddef>
#include <numeric>
#include "Group.h"

// Projectors onto one-dimensional representations are uniquely specified using 
// two characters: the rotation character (+1 or -1) and the reflection character (+1 or -1). 
std::vector<std::vector<double>> one_dimension_projector_generator(
    size_t order, double rot_char, double ref_char) {
    std::vector<std::vector<double>> projectors;
    projectors.resize(1, std::vector<double>(2 * order, 0.0));
    projectors[0][0] = 1;
    for (int rot = 1; rot < order; ++rot) {
        projectors[0][rot] = std::pow(rot_char, rot);
    }
    for (int ref = 0; ref < order; ++ref) {
        projectors[0][order + ref] = ref_char * std::pow(rot_char, ref);
    }
    return std::move(projectors);
}

// Coefficients of projectors onto two-dimensional representations are elements
// of two-by-two matrices of rotation and reflactions. 
std::vector<std::vector<double>> two_dimension_projector_generator(
    size_t order, size_t two_dim_repr) {
    std::vector<std::vector<double>> projectors;
    projectors.resize(4, std::vector<double>(2 * order, 0.0));
    // identity:
    projectors[0][0] = 1.0;
    projectors[1][0] = 0.0;
    projectors[2][0] = 0.0;
    projectors[3][0] = 1.0;
    // rotations:
    for (int rot = 1; rot < order; ++rot) {
        double angle = 2 * M_PI * (two_dim_repr + 1) * rot / order;
        projectors[0][rot] = std::cos(angle);
        projectors[1][rot] = std::sin(angle);
        projectors[2][rot] = -1.0 * std::sin(angle);
        projectors[3][rot] = std::cos(angle);
    }
    // reflections:
    for (int ref = 0; ref < order; ++ref) {
        double angle = 2 * M_PI * (two_dim_repr + 1) * ref / order;
        projectors[0][order + ref] = std::cos(angle);
        projectors[1][order + ref] = std::sin(angle);
        projectors[2][order + ref] = std::sin(angle);
        projectors[3][order + ref] = -1.0 * std::cos(angle);
    }
    for (auto& projector : projectors) {
        for (auto& el : projector) {
            if (std::abs(el) < 1e-6) {
                el = 0.0;
            }
        }
    }
    return std::move(projectors);
}

group::Group::AlgebraicProperties constructDihedral(unsigned int order)
{
    if (order == 0) {
        throw group::InitializationError("Order cannot be equal to zero");
    }
    if (order == 1 || order == 2) {
        throw group::InitializationError("For construction of Dih1 and Dih2 groups use S2 group");
    }
    group::Group::AlgebraicProperties properties;
    properties.group_size = 2 * order;
    properties.is_abelian = false;
    // generators are one rotation and one mirror:
    properties.number_of_generators = 2;
    properties.orders_of_generators = {order, 2};

    for (size_t ref = 0; ref < 2; ++ref) {
        for (size_t rot = 0; rot < order; ++rot) {
            properties.group_in_form_of_generators.emplace_back(std::vector<size_t>({rot, ref}));
        }
    }

    // identity:
    properties.orders_of_elements.push_back(1);
    // rotations:
    for (int rot = 1; rot < order; ++rot) {
        properties.orders_of_elements.push_back(order / std::gcd(rot, order));
    }
    // reflections:
    for (int ref = order; ref < 2 * order; ++ref) {
        properties.orders_of_elements.push_back(2);
    }

    if (order % 2 == 0) {
        properties.number_of_representations = (order + 6) / 2;
        properties.coefficients_of_projectors.resize(properties.number_of_representations);
        // four one-dimension representation:
        for (int repr = 0; repr < 4; ++repr) {
            properties.dimension_of_representation.emplace_back(1);
            properties.number_of_projectors_of_representation.emplace_back(1);
        }
        // (order/2 - 1) two-dimension representation:
        for (int repr = 4; repr < (order + 6) / 2; ++repr) {
            properties.dimension_of_representation.emplace_back(2);
            properties.number_of_projectors_of_representation.emplace_back(4);
            properties.coefficients_of_projectors[repr].resize(4);
        }
        // A1:
        properties.coefficients_of_projectors[0] = 
            one_dimension_projector_generator(order, 1.0, 1.0);
        // A2:
        properties.coefficients_of_projectors[1] = 
            one_dimension_projector_generator(order, 1.0, -1.0);;
        // B1
        properties.coefficients_of_projectors[2] = 
            one_dimension_projector_generator(order, -1.0, 1.0);
        // B2
        properties.coefficients_of_projectors[3] = 
            one_dimension_projector_generator(order, -1.0, -1.0);
        // all Es:
        for (size_t two_dim_repr = 0; two_dim_repr < (order/2 - 1); ++two_dim_repr) {
            properties.coefficients_of_projectors[4 + two_dim_repr] = 
                two_dimension_projector_generator(order, two_dim_repr);
        }
    } else {
        properties.number_of_representations = (order + 3) / 2;
        properties.coefficients_of_projectors.resize(properties.number_of_representations);
        // two one-dimension representation:
        for (int repr = 0; repr < 2; ++repr) {
            properties.dimension_of_representation.emplace_back(1);
            properties.number_of_projectors_of_representation.emplace_back(1);
        }
        // (order - 1)/2 two-dimension representation:
        for (int repr = 2; repr < (order + 3) / 2; ++repr) {
            properties.dimension_of_representation.emplace_back(2);
            properties.number_of_projectors_of_representation.emplace_back(4);
        }
        // A1:
        properties.coefficients_of_projectors[0] = 
            one_dimension_projector_generator(order, 1.0, 1.0);
        // A2:
        properties.coefficients_of_projectors[1] = 
            one_dimension_projector_generator(order, 1.0, -1.0);;
        // all Es:
        for (size_t two_dim_repr = 0; two_dim_repr < (order - 1)/2; ++two_dim_repr) {
            properties.coefficients_of_projectors[2 + two_dim_repr] = 
                two_dimension_projector_generator(order, two_dim_repr);
        }
    }

    // TODO: properties.cayley_table;
    return properties;
}


// clang-format off

/*
 It is symmetric group of order two. It has two group elements, two representations.
 Maximum size of orbit -- two, so the group has two projectors, one for each representation.
 */
const group::Group::AlgebraicProperties GroupInfoS2 = {2,
                                                       2,
                                                       true,
                                                       {1, 1},
                                                       {1, 1},
                                                       1,
                                                       {2},
                                                       {{0}, {1}},
                                                       {1, 2},
                                                       {{{1,  1}},  // "a" representation
                                                       {{1, -1}}},  // "b" representation
                                                       {{{0, 0}, {0}}, // "a" * "a" = "a"
                                                       {{0, 1}, {1}},  // "a" * "b" = "b"
                                                       {{1, 0}, {1}},  // "b" * "a" = "b"
                                                       {{1, 1}, {0}}}  // "b" * "b" = "a"
};

/*
 It is symmetric group of order three. It has six group elements, three representations: two one-dimensional
 and one two-dimensional.
 Maximum size of orbit -- six, so the group has six projectors, one for each one-dimensional representation
 and four for two-dimensional representations.
 */
const group::Group::AlgebraicProperties GroupInfoS3 = {6,
                                                       3,
                                                       false,
                                                       {1, 1, 2},
                                                       {1, 1, 4},
                                                       2,
                                                       {3, 2},
                                                       {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}},
                                                       {1, 3, 3, 2, 2, 2},
                                                       {{{1,  1,  1,  1,  1,  1}}, // "a" representation
                                                        {{1,  1,  1, -1, -1, -1}}, // "b" representation
                                                        {{2, -1, -1,  2,  -1,  -1},
                                                        {2,  -1,  -1,  -2, 1, 1},
                                                        {0,  1,  -1,  0,  -1, 1},
                                                        {0,  1, -1,  0,  1,  -1}}},  // "e" representation
                                                        {{{0, 0}, {0}},      // "a" * "a" = "a"
                                                        {{0, 1}, {1}},       // "a" * "b" = "b"
                                                        {{0, 2}, {2}},       // "a" * "e" = "e"
                                                        {{1, 0}, {1}},       // "b" * "a" = "b"
                                                        {{1, 1}, {0}},       // "b" * "b" = "a"
                                                        {{1, 2}, {2}},       // "b" * "e" = "e"
                                                        {{2, 0}, {2}},       // "e" * "a" = "e"
                                                        {{2, 1}, {2}},       // "e" * "b" = "e"
                                                        {{2, 2}, {0, 1, 2}}} // "e" * "e" = "a" + "b" + "e"
};

/*
 It is subgroup of symmetric group of order four. It has eight group elements, five representations: four one-dimensional
 and one two-dimensional. This group is isomorphic to D4 point group and represents four-square graph.
 Maximum size of orbit -- eight, so the group has eight projectors, one for each one-dimensional representation
 and four for two-dimensional representations.
 */
 const group::Group::AlgebraicProperties GroupInfoD4 = {8,
                                                        5,
                                                        false,
                                                        {1, 1, 1, 1, 2},
                                                        {1, 1, 1, 1, 4},
                                                        2,
                                                        {4, 2},
                                                        {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {0, 1}, {1, 1}, {2, 1}, {3, 1}},
                                                        {1, 4, 2, 4, 2, 2, 2, 2},
                                                        {{{1,  1,  1,  1,  1,  1,  1,  1}}, // "a1" representation
                                                        {{1,  1,  1,  1, -1, -1, -1, -1}}, // "a2" representation
                                                        {{1, -1,  1, -1,  1, -1,  1, -1}}, // "b1" representation
                                                        {{1, -1,  1, -1, -1,  1, -1,  1}}, // "b2" representation
                                                        {{2,2,-2,-2,2,2,-2,-2},
                                                        {2,-2,-2,2,2,-2,-2,2},
                                                        {2,-2,-2,2,-2,2,2,-2},
                                                        {2,2,-2,-2,-2,-2,2,2}}},  // "e" representation
                                                        {{{0, 0}, {0}},      // "a1" * "a1" = "a1"
                                                        {{0, 1}, {1}},       // "a1" * "a2" = "a2"
                                                        {{0, 2}, {2}},       // "a1" * "b1" = "b1"
                                                        {{0, 3}, {3}},       // "a1" * "b2" = "b2"
                                                        {{0, 4}, {4}},       // "a1" * "e" = "e"
                                                        {{1, 0}, {1}},       // "a2" * "a1" = "a2"
                                                        {{1, 1}, {0}},       // "a2" * "a2" = "a1"
                                                        {{1, 2}, {3}},       // "a2" * "b1" = "b2"
                                                        {{1, 3}, {2}},       // "a2" * "b2" = "b1"
                                                        {{1, 4}, {4}},       // "a2" * "e" = "e"
                                                        {{2, 0}, {2}},       // "b1" * "a1" = "b1"
                                                        {{2, 1}, {3}},       // "b1" * "a2" = "b2"
                                                        {{2, 2}, {0}},       // "b1" * "b1" = "a1"
                                                        {{2, 3}, {1}},       // "b1" * "b2" = "a2"
                                                        {{2, 4}, {4}},       // "b1" * "e" = "e"
                                                        {{3, 0}, {3}},       // "b2" * "a1" = "b2"
                                                        {{3, 1}, {2}},       // "b2" * "a2" = "b1"
                                                        {{3, 2}, {1}},       // "b2" * "b1" = "a2"
                                                        {{3, 3}, {0}},       // "b2" * "b2" = "a1"
                                                        {{3, 4}, {4}},       // "b2" * "e" = "e"
                                                        {{4, 0}, {4}},       // "e" * "a1" = "e"
                                                        {{4, 1}, {4}},       // "e" * "a2" = "e"
                                                        {{4, 2}, {4}},       // "e" * "b1" = "e"
                                                        {{4, 3}, {4}},       // "e" * "b2" = "e"
                                                        {{4, 4}, {0, 1, 2, 3}},       // "e" * "e" = "a1" + "a2" + "b1" + "b2"
}
};

// clang-format on