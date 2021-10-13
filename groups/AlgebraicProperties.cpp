#include "Group.h"

/*
 It is symmetric group of order two. It has two group elements, two representations.
 Maximum size of orbit -- two, so the group has two projectors, one for each representation.
 */
const Group::AlgebraicProperties GroupInfoS2 = {2,
                                                2,
                                                true,
                                                {1, 1},
                                                {1, 1},
                                                1,
                                                {2},
                                                {{0}, {1}},
                                                {1, 2},
                                                {{{1,  1}},  // "a" representation
                                                {{1, -1}}}  // "b" representation
};

/*
 It is symmetric group of order three. It has six group elements, three representations: two one-dimensional
 and one two-dimensional.
 Maximum size of orbit -- six, so the group has six projectors, one for each one-dimensional representation
 and four for two-dimensional representations.
 */
const Group::AlgebraicProperties GroupInfoS3 = {6,
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
                                                 {0,  0,  0,  0,  1, -1}}}  // "e" representation
};