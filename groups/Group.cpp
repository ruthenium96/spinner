#include "Group.h"

//Group::Group(GroupNames group_name, std::vector<Permutation> generators) :
//generators_(std::move(generators)), groupInfo(return_group_info_by_group_name(group_name)) {
//    Permutation identity(generators_[0].size());
//    for (uint32_t i = 0; i < generators_[0].size(); ++i) {
//        identity[i] = i;
//    }
//    elements_.resize(groupInfo.group_size);
//    for (uint32_t i = 0; i < groupInfo.group_size; ++i) {
//        Permutation element = identity;
//        // по элементам группы
//        for (uint32_t j = 0; j < groupInfo.number_of_generators; ++j) {
//            // по генераторам группы
//            for (uint32_t k = 0; k < groupInfo.group_in_form_of_generators[i][j]; ++k) {
//                // по степени генератора
//                Permutation result_element(element.size());
//                for (uint32_t l = 0; l < generators_[0].size(); ++l) {
//                    // по передвигаемым спинам
//                    result_element[l] = element[generators_[j][l]];
//                }
//                element = result_element;
//            }
//        }
//        elements_[i] = element;
//    }
//}

std::vector<std::vector<Projection>> Group::permutate(const std::vector<Projection> &initial) const {
    std::vector<std::vector<Projection>> permutated_vectors(groupInfo.group_size);
    for (size_t i = 0; i < groupInfo.group_size; ++i) {
        permutated_vectors[i].resize(initial.size());
        for (size_t j = 0; j < initial.size(); ++j) {
            size_t position = elements_[i][j];
            permutated_vectors[i][position] = initial[j];
        }
    }
    return std::move(permutated_vectors);
}

Group::Group(const GroupInfo& group_info, std::vector<Permutation> generators) :
generators_(std::move(generators)), groupInfo(group_info) {
    Permutation identity(generators_[0].size());
    for (uint32_t i = 0; i < generators_[0].size(); ++i) {
        identity[i] = i;
    }
    elements_.resize(groupInfo.group_size);
    // cycle over group elements
    for (uint32_t i = 0; i < groupInfo.group_size; ++i) {
        Permutation element = identity;
        // cycle over different type of generators
        for (uint32_t j = 0; j < groupInfo.number_of_generators; ++j) {
            // cycle over degree of generator
            for (uint32_t k = 0; k < groupInfo.group_in_form_of_generators[i][j]; ++k) {
                Permutation result_element(element.size());
                // do the permutation:
                for (uint32_t l = 0; l < generators_[0].size(); ++l) {
                    result_element[l] = element[generators_[j][l]];
                }
                element = result_element;
            }
        }
        elements_[i] = element;
    }
}
