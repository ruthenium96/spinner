#ifndef JULY_GROUP_INFO_H
#define JULY_GROUP_INFO_H

// TODO: create simple way to keep _all_ projectors for non-Abelian groups

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
    std::vector<std::vector<double>> coefficients_of_projectors;
};

static const GroupInfo GroupInfoP2 = {2, // group size
                                      2, // number of representations
                                      1, // number of generators
                                      // group elements as product of generators:
                                      {{0}, {1}},
                                      // coefficients of projectors
                                      {{1, 1},  // "a" representation
                                       {1, -1}} // "b" representation
};

static const GroupInfo GroupInfoP3 = {6, // group size
                                      3, // number of representations
                                      2, // number of generators
                                      // group elements as product of generators:
                                      {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {2, 1}},
                                      // coefficients of projectors
                                      {{1, 1, 1, 1, 1, 1},    // "a" representation
                                       {1, 1, 1, -1, -1, -1}, // "b" representation
                                       {2, -1, -1, 0, 0, 0}}  // "e" representation
};

enum GroupNames {
    P2,
    P3,
};

static const GroupInfo& return_group_info_by_group_name(GroupNames group_name) {
    if (group_name == P2) {
        return GroupInfoP2;
    }
    if (group_name == P3) {
        return GroupInfoP3;
    }
}


#endif //JULY_GROUP_INFO_H
