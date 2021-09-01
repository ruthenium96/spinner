#ifndef JULY_GROUP_INFO_H
#define JULY_GROUP_INFO_H

struct GroupInfo {
    // the number of permutations in the group
    uint8_t group_size;
    // the number of representations, also the number of conjugacy class
    uint8_t number_of_representations;
    // the minimal number of elements required to produce all elements of the group
    uint8_t number_of_generators;
    // all elements of the group in the form a^vector[0] * b^vector[1] * ...,
    // where {a, b, ...} -- generators in the _descending_ order of element order.
    // element order is the minimum m, such that g^m = e
    // It has the size (group_size) x (number of generators).
    std::vector<std::vector<size_t>> group_in_form_of_generators;
    // the coefficients of projectors to representations
    // It has the size (number_of_representations) x (group_size).
    std::vector<std::vector<Coefficient>> coefficients_of_projectors;
};

const GroupInfo GroupInfoP2 = {2, 2, 1,
                               {{0}, {1}},
                               {{1, 1},
                                {1, -1}}
};

const GroupInfo GroupInfoP3 = {6, 3, 2,
                               {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}},
                               {{1, 1, 1, 1, 1, 1},
                                {1, 1, 1, -1, -1, -1},
                                {2, -1, -1, 0, 0, 0}}
};

enum GroupNames {
    P2,
    P3,
};

//const GroupInfo& return_group_info_by_group_name(GroupNames group_name) {
//    if (group_name == P2) {
//        return GroupInfoP2;
//    }
//    if (group_name == P3) {
//        return GroupInfoP3;
//    }
//}


#endif //JULY_GROUP_INFO_H
