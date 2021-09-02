#include "Group.h"

Group::Group(GroupNames group_name, std::vector<Permutation> generators) :
        generators_(std::move(generators)), info(return_group_info_by_group_name(group_name)) {
    Permutation identity(generators_[0].size());
    for (uint32_t i = 0; i < generators_[0].size(); ++i) {
        identity[i] = i;
    }
    elements_.resize(info.group_size);
    // over group elements
    for (uint32_t i = 0; i < info.group_size; ++i) {
        Permutation element = identity;
        // over group generators
        for (uint32_t j = 0; j < info.number_of_generators; ++j) {
            // over generator's power
            for (uint32_t k = 0; k < info.group_in_form_of_generators[i][j]; ++k) {
                Permutation result_element(element.size());
                // over spin centers
                for (uint32_t l = 0; l < generators_[0].size(); ++l) {
                    result_element[l] = element[generators_[j][l]];
                }
                element = result_element;
            }
        }
        elements_[i] = element;
    }
}

std::vector<std::vector<uint8_t>> Group::permutate(const std::vector<uint8_t> &initial) const {
    std::vector<std::vector<uint8_t>> permutated_vectors(info.group_size);
    for (size_t i = 0; i < info.group_size; ++i) {
        permutated_vectors[i].resize(initial.size());
        for (size_t j = 0; j < initial.size(); ++j) {
            size_t position = elements_[i][j];
            permutated_vectors[i][position] = initial[j];
        }
    }
    return std::move(permutated_vectors);
}