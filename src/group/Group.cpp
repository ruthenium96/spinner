#include "Group.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace {

bool NumberOfGeneratorsConsistent(
    std::vector<group::Permutation> const& generators,
    group::Group::AlgebraicProperties const& info) {
    return generators.size() == info.number_of_generators;
}

bool SizesAreEqual(
    std::vector<group::Permutation> const& generators,
    group::Group::AlgebraicProperties const& info) {
    for (size_t i = 0; i < info.number_of_generators; ++i) {
        if (generators[0].size() != generators[i].size()) {
            return false;
        }
    }
    return true;
}

bool GeneratorIsValid(
    std::vector<group::Permutation> const& generators,
    group::Group::AlgebraicProperties const& info) {
    for (size_t i = 0; i < info.number_of_generators; ++i) {
        group::Permutation generator_sorted = generators[i];
        std::sort(generator_sorted.begin(), generator_sorted.end());
        for (size_t j = 0; j < generator_sorted.size(); ++j) {
            if (generator_sorted[j] != j) {
                return false;
            }
        }
    }
    return true;
}

group::Permutation ElementInPower(const group::Permutation& element, size_t power) {
    group::Permutation result = element;
    for (size_t i = 1; i < power; ++i) {
        group::Permutation tmp_result(result.size());
        for (size_t n = 0; n < tmp_result.size(); ++n) {
            tmp_result[n] = result[element[n]];
        }
        result = tmp_result;
    }
    return result;
}

bool ElementsInPowerOfItsOrderIsIdentity(
    const std::vector<group::Permutation>& elements,
    const std::vector<size_t>& order_of_elements) {
    group::Permutation identity(elements[0].size());
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

}  // namespace

namespace group {
Group::Group(GroupTypeEnum group_name, std::vector<Permutation> generators) :
    generators_(std::move(generators)),
    properties(Group::return_group_info_by_group_name(group_name)) {
    if (!NumberOfGeneratorsConsistent(generators_, properties)) {
        throw InitializationError(
            "The number of generators does not equal to the number of group number_of_generators.");
    }

    if (!SizesAreEqual(generators_, properties)) {
        throw InitializationError("The sizes of generators are different.");
    }

    if (!GeneratorIsValid(generators_, properties)) {
        throw InitializationError("Generator is invalid");
    }

    if (!ElementsInPowerOfItsOrderIsIdentity(generators_, properties.orders_of_generators)) {
        throw InitializationError("Generator in power of its order does not equal identity.");
    }

    Permutation identity(generators_[0].size());
    for (size_t i = 0; i < generators_[0].size(); ++i) {
        identity[i] = i;
    }
    elements_.resize(properties.group_size);
    // over group elements
    for (uint32_t i = 0; i < properties.group_size; ++i) {
        Permutation element = identity;
        // over group generators
        for (uint32_t j = 0; j < properties.number_of_generators; ++j) {
            // over generator's power
            for (uint32_t k = 0; k < properties.group_in_form_of_generators[i][j]; ++k) {
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

    if (!ElementsInPowerOfItsOrderIsIdentity(elements_, properties.orders_of_elements)) {
        throw InitializationError("Element in power of its order does not equal identity.");
    }
}

std::vector<std::vector<uint8_t>> Group::permutate(const std::vector<uint8_t>& initial) const {
    std::vector<std::vector<uint8_t>> permutated_vectors(properties.group_size);
    for (size_t i = 0; i < properties.group_size; ++i) {
        permutated_vectors[i].resize(initial.size());
        for (size_t j = 0; j < initial.size(); ++j) {
            size_t position = elements_[i][j];
            permutated_vectors[i][position] = initial[j];
        }
    }
    return permutated_vectors;
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

const Group::AlgebraicProperties&
Group::return_group_info_by_group_name(Group::GroupTypeEnum group_name) {
    if (group_name == S2) {
        return GroupInfoS2;
    }
    if (group_name == S3) {
        return GroupInfoS3;
    }
}

const std::vector<Permutation>& Group::getElements() const {
    return elements_;
}
}  // namespace group