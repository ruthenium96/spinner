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

bool ElementsInPowerLowerThanItsOrderIsNotIdentity(
    const std::vector<group::Permutation>& elements,
    const std::vector<size_t>& order_of_elements) {
    group::Permutation identity(elements[0].size());
    for (size_t i = 0; i < elements[0].size(); ++i) {
        identity[i] = i;
    }
    for (size_t g = 0; g < elements.size(); ++g) {
        for (int order = 1; order < order_of_elements[g]-1; ++order) {
            if (ElementInPower(elements[g], order) == identity) {
                return false;
            }    
        }
    }
    return true;
}

}  // namespace

namespace group {
Group::Group(GroupType group_type, std::vector<Permutation> generators) :
    generators_(std::move(generators)),
    properties(Group::return_group_info_by_group_type(group_type)) {
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

    if (!ElementsInPowerLowerThanItsOrderIsNotIdentity(generators_, properties.orders_of_generators)) {
        throw InitializationError("Generator in power lower than its order equals to identity.");
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
    
    if (!ElementsInPowerLowerThanItsOrderIsNotIdentity(elements_, properties.orders_of_elements)) {
        throw InitializationError("Element in power lower than its order equals to identity.");
    }
}

// Returns the indexes of multiplicity centers, forming the same orbit.
std::vector<std::set<size_t>> Group::construct_orbits_of_mults() const {
    std::vector<std::set<size_t>> orbits_of_mults;

    size_t number_of_mults = elements_[0].size();

    std::vector<bool> does_mult_already_formed_orbit(number_of_mults, false);

    for (size_t i = 0; i < number_of_mults; ++i) {
        std::set<size_t> current_orbit;
        for (size_t j = 0; j < properties.group_size; ++j) {
            size_t current_position = elements_[j][i];
            if (!does_mult_already_formed_orbit[current_position]) {
                current_orbit.insert(current_position);
                does_mult_already_formed_orbit[current_position] = true;
            }
        }
        if (!current_orbit.empty()) {
            orbits_of_mults.emplace_back(std::move(current_orbit));
        }

        if (std::all_of(
                does_mult_already_formed_orbit.begin(),
                does_mult_already_formed_orbit.end(),
                [](bool b) { return b; })) {
            break;
        }
    }

    return orbits_of_mults;
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

Group::AlgebraicProperties
Group::return_group_info_by_group_type(Group::GroupType group_type) {
    if (group_type.type_enum == S2) {
        return GroupInfoS2;
    }
    if (group_type.type_enum == Dihedral) {
        if (!group_type.order.has_value()) {
            throw InitializationError("Dihedral must have an order");
        }
        return constructDihedral(group_type.order.value());
    }
    throw InitializationError("Unknown group_type" + std::to_string(group_type.type_enum));
}

const std::vector<Permutation>& Group::getElements() const {
    return elements_;
}

const std::vector<Permutation>& Group::getGenerators() const {
    return generators_;
}


bool Group::do_groups_commute(const Group& rhs) const {
    size_t lhs_size = generators_.size();
    size_t rhs_size = rhs.generators_.size();
    size_t permutation_size = elements_.back().size();

    for (size_t l = 0; l < lhs_size; ++l) {
        const auto& lg = generators_.at(l);
        for (size_t r = 0; r < rhs_size; ++r) {
            const auto& rg = rhs.generators_.at(r);
            for (size_t el = 0; el < permutation_size; ++el) {
                if (rg.at(lg.at(el)) != lg.at(rg.at(el))) {
                    return false;
                }
            }
        }
    }
    return true;
}

size_t Group::size_of_permutations() const {
    return elements_.back().size();
}
}  // namespace group