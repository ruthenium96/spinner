#include "Group.h"
#include <stdexcept>
#include <algorithm>
#include <numeric>

Group::Group(group::GroupNames group_name, std::vector<Permutation> generators) :
        generators_(std::move(generators)), info(group::return_group_info_by_group_name(group_name)) {
    if (generators_.size() != info.number_of_generators) {
        throw std::length_error("The number of generators does not equal to the number of group number_of_generators.");
    }
    for (size_t i = 0; i < info.number_of_generators; ++i) {
        if (generators_[0].size() != generators_[i].size()) {
            throw std::length_error("The sizes of generators are different.");
        }
        std::vector<bool> contains_number(generators_[i].size(), false);
        for (unsigned char j : generators_[i]) {
            if (contains_number[j]) {
                throw std::invalid_argument("Generator contains number twice.");
            } else {
                contains_number[j] = true;
            }
        }
        if (std::find(contains_number.begin(), contains_number.end(), false) != contains_number.end()) {
            throw std::invalid_argument("Generator does not contain all numbers.");
        }
    }
    // TODO: check that generator ^ {order_of_generator} == identity
    //  how to find / save orders of generators?


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

    // TODO: check that element ^ {order_of_element} == identity
    //  how to find / save orders of elements?

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