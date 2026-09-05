#include "OrderOfSummation.h"

#include <algorithm>
#include <optional>
#include <set>
#include <stdexcept>


namespace {
size_t find_root(size_t pos, const std::vector<std::optional<size_t>>& parent_of_node) {
    while (parent_of_node.at(pos).has_value()) {
        pos = parent_of_node.at(pos).value();
    }
    return pos;
}
} // namespace

namespace spinner::index_converter::s_squared {

std::shared_ptr<const OrderOfSummation> OrderOfSummation::constructFromOrbits(
    const std::vector<std::vector<std::set<size_t>>>& all_groups_orbits_of_mults,
    size_t number_of_mults) {
    if (number_of_mults == 0) {
        throw std::invalid_argument("Cannot construct order of summation for zero spins");
    }

    const size_t number_of_summation = number_of_mults - 1;
    const size_t number_of_nodes = number_of_mults + number_of_summation;

    size_t performed_summations = 0;

    // For every node of the coupling tree, stores its parent.
    // std::nullopt means that the node is currently a root.
    std::vector<std::optional<size_t>> parent_of_node(number_of_nodes, std::nullopt);

    // Roots of the currently disconnected coupling subtrees.
    std::set<size_t> active_roots;
    for (size_t i = 0; i < number_of_mults; ++i) {
        active_roots.insert(i);
    }

    auto order_of_summations = std::make_shared<OrderOfSummation>();

    auto add_summation =
        [&](size_t left, size_t right) {
            if (left == right) {
                throw std::logic_error("Cannot sum a coupling subtree with itself");
            }

            if (active_roots.find(left) == active_roots.end()
                || active_roots.find(right) == active_roots.end()) {
                throw std::logic_error("Attempt to sum a coupling subtree which is not active");
            }

            const size_t result = number_of_mults + performed_summations;

            AdditionInstruction instruction;
            instruction.left = left;
            instruction.right = right;
            instruction.result = result;

            parent_of_node.at(left) = result;
            parent_of_node.at(right) = result;

            active_roots.erase(left);
            active_roots.erase(right);
            active_roots.insert(result);

            order_of_summations->instructions_.push_back(instruction);

            ++performed_summations;
        };

    for (size_t number_of_group = 0;
         number_of_group < all_groups_orbits_of_mults.size() && performed_summations < number_of_summation;
         ++number_of_group) {
        const auto& current_orbits = all_groups_orbits_of_mults.at(number_of_group);

        for (const auto& orbit : current_orbits) {
            if (performed_summations == number_of_summation) {
                break;
            }

            if (orbit.size() == 1) {
                continue;
            }

            std::vector<size_t> positions_of_summands;

            for (const auto pos : orbit) {
                const size_t root = find_root(pos, parent_of_node);

                if (std::find(
                        positions_of_summands.begin(),
                        positions_of_summands.end(),
                        root)
                    == positions_of_summands.end()) {
                    positions_of_summands.push_back(root);
                }
            }

            if (positions_of_summands.size() == 1) {
                // All elements of the orbit already belong
                // to the same coupling subtree.
                continue;
            }

            if (positions_of_summands.size() != 2) {
                throw std::invalid_argument(
                    "Cannot construct binary order of summation "
                    "from an orbit with more than two active subspaces");
            }

            add_summation(positions_of_summands[0], positions_of_summands[1]);
        }
    }

    // If the symmetry information has not determined the complete
    // coupling tree, join the remaining subtrees pairwise.
    while (performed_summations < number_of_summation) {
        if (active_roots.size() < 2) {
            throw std::logic_error("Cannot find two active subspaces for the next summation");
        }

        auto root_it = active_roots.begin();
        const size_t left = *root_it;
        ++root_it;
        const size_t right = *root_it;

        add_summation(left, right);
    }

    if (performed_summations != number_of_summation) {
        throw std::logic_error("Incorrect number of summations in the constructed coupling tree");
    }

    if (active_roots.size() != 1) {
        throw std::logic_error("Constructed order of summation does not form a single coupling tree");
    }

    const size_t expected_root = number_of_nodes - 1;

    if (*active_roots.begin() != expected_root) {
        throw std::logic_error("Unexpected root index in the constructed coupling tree");
    }

    return order_of_summations;
}

std::vector<OrderOfSummation::AdditionInstruction>::const_iterator OrderOfSummation::begin() const {
    return instructions_.begin();
}

std::vector<OrderOfSummation::AdditionInstruction>::const_iterator OrderOfSummation::end() const {
    return instructions_.end();
}

const OrderOfSummation::AdditionInstruction& OrderOfSummation::at(size_t i) const {
    return instructions_.at(i);
}

size_t OrderOfSummation::size() const {
    return instructions_.size();
}

const std::vector<OrderOfSummation::AdditionInstruction>& OrderOfSummation::getInstructions() const {
    return instructions_;
}


} // namespace spinner::index_converter::s_squared