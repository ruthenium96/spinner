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
                                                        {{2, -1, -1,  0,  0,  0},
                                                         {0,  1, -1,  0,  0,  0},
                                                         {0,  0,  0,  2, -1, -1},
                                                         {0,  0,  0,  0,  1, -1}}},  // "e" representation
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

// clang-format on