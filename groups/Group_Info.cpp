#include "Group_Info.h"

namespace group {

    /*
     It is symmetric group of order two. It has two group elements, two representations.
     Maximum size of orbit -- two, so the group has two projectors, one for each representation.
     */
    const GroupInfo GroupInfoS2 = {2,
                                   2,
                                   {1, 1},
                                   {1, 1},
                                   1,
                                   {{0}, {1}},
                                   {{{1,  1}},  // "a" representation
                                    {{1, -1}}}  // "b" representation
    };

    /*
     It is symmetric group of order three. It has six group elements, three representations: two one-dimensional
     and one two-dimensional.
     Maximum size of orbit -- six, so the group has six projectors, one for each one-dimensional representation
     and four for two-dimensional representations.
     */
    const GroupInfo GroupInfoS3 = {6,
                                   3,
                                   {1, 1, 2},
                                   {1, 1, 4},
                                   2,
                                   {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}},
                                   {{{1,  1,  1,  1,  1,  1}}, // "a" representation
                                    {{1,  1,  1, -1, -1, -1}}, // "b" representation
                                    {{2, -1, -1,  0,  0,  0},
                                     {0,  1, -1,  0,  0,  0},
                                     {0,  0,  0,  2, -1, -1},
                                     {0,  0,  0,  0,  1, -1}}}  // "e" representation
    };

    const GroupInfo& return_group_info_by_group_name(GroupNames group_name) {
        if (group_name == S2) {
            return GroupInfoS2;
        }
        if (group_name == S3) {
            return GroupInfoS3;
        }
    }
}