#include "Group.h"
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
                                                        {{2,  0, -2,  0,  0,  0,  0,  0},
                                                        {0,  2,  0, -2,  0,  0,  0,  0},
                                                        {0,  0,  0,  0,  2,  0, -2,  0},
                                                        {0,  0,  0,  0,  0,  2,  0, -2}}},  // "e" representation
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