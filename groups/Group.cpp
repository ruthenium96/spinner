#include "Group.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace {

bool NumberOfGeneratorsConsistent(std::vector<Permutation> const& generators,
                             group::GroupInfo const& info) {
    return generators.size() == info.number_of_generators;
}

bool SizesAreEqual(std::vector<Permutation> const& generators, group::GroupInfo const& info) {
    for (size_t i = 0; i < info.number_of_generators; ++i) {
        if (generators[0].size() != generators[i].size()) {
            return false;
        }
    }
    return true;
}

bool GeneratorIsValid(std::vector<Permutation> const& generators, group::GroupInfo const& info) {
    for (size_t i = 0; i < info.number_of_generators; ++i) {
        Permutation generator_sorted = generators[i];
        std::sort(generator_sorted.begin(), generator_sorted.end());
        for (size_t j = 0; j < generator_sorted.size(); ++j) {
            if (generator_sorted[j] != j) {
                return false;
            }
        }
    }
    return true;
}

Permutation ElementInPower(const Permutation& element, size_t power) {
    Permutation result = element;
    for (size_t i = 1; i < power; ++i) {
        Permutation tmp_result(result.size());
        for (size_t n = 0; n < tmp_result.size(); ++n) {
            tmp_result[n] = result[element[n]];
        }
        result = tmp_result;
    }
    return std::move(result);
};

bool ElementsInPowerOfItsOrderIsIdentity(const std::vector<Permutation>& elements,
                                         const std::vector<size_t>& order_of_elements) {
    Permutation identity(elements[0].size());
    for (size_t i = 0; i < elements[0].size(); ++i) {
        identity[i] = i;
    }
    for (size_t i = 0; i < elements.size(); ++i) {
        if (ElementInPower(elements[i], order_of_elements[i]) != identity) {
            return false;
        }
    }
    return true;
}

} // namespace

Group::Group(group::GroupNames group_name, std::vector<Permutation> generators)
    : generators_(std::move(generators)), info(group::return_group_info_by_group_name(group_name)) {

    if (!NumberOfGeneratorsConsistent(generators_, info)) {
        throw std::length_error(
            "The number of generators does not equal to the number of group number_of_generators.");
    }

    if (!SizesAreEqual(generators_, info)) {
        throw std::length_error("The sizes of generators are different.");
    }

    if (!GeneratorIsValid(generators_, info)) {
        throw std::invalid_argument("Generator is invalid.");
    }

    if (!ElementsInPowerOfItsOrderIsIdentity(generators_, info.orders_of_generators)) {
        throw std::invalid_argument("Generator in power of its order does not equal identity.");
    }

    Permutation identity(generators_[0].size());
    for (size_t i = 0; i < generators_[0].size(); ++i) {
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

    if (!ElementsInPowerOfItsOrderIsIdentity(elements_, info.orders_of_elements)) {
        throw std::invalid_argument("Element in power of its order does not equal identity.");
    }
}

std::vector<std::vector<uint8_t>> Group::permutate(const std::vector<uint8_t>& initial) const {
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

/*
 Different sets of generators can produce isomorphic group.
 So we cannot compare the generators.
 The isomorphic groups in our case differ only in elements order,
 thus sorting will help to unify them.
 */
bool Group::operator==(const Group& rhs) const {
    std::vector<Permutation> left_elements_copy = elements_;
    std::vector<Permutation> right_elements_copy = rhs.elements_;
    std::sort(left_elements_copy.begin(), left_elements_copy.end());
    std::sort(right_elements_copy.begin(), right_elements_copy.end());
    return left_elements_copy == right_elements_copy;
}

bool Group::operator!=(const Group& rhs) const {
    return !(rhs == *this);
}
