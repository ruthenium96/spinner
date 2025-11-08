#include "OptimizationsParser.h"
#include <cstddef>
#include <optional>
#include "Tools.h"
#include "src/common/physical_optimization/OptimizationList.h"

namespace input {
OptimizationsParser::OptimizationsParser(YAML::Node optimizations_node) {
    auto mode_string = extractValue<std::string>(optimizations_node, "mode");
    if (mode_string == "none") {
        optimizations_list_ = common::physical_optimization::OptimizationList();
        // do nothing
    } else if (mode_string == "custom") {
        customParser(extractValue<YAML::Node>(optimizations_node, "custom"));
    } else if (mode_string == "auto") {
        throw std::invalid_argument("optimizations::mode == auto is not implemented yet");
    } else {
        throw std::invalid_argument("Incorrect argument of optimizations::mode: " + mode_string);
    }
    throw_if_node_is_not_empty(optimizations_node);
}

const std::optional<common::physical_optimization::OptimizationList>&
OptimizationsParser::getOptimizationList() const {
    return optimizations_list_;
}

void OptimizationsParser::customParser(YAML::Node custom_node) {
    auto basis_type_ = extractValue<common::physical_optimization::OptimizationList::BasisType>(custom_node, "basis");
    optimizations_list_ = common::physical_optimization::OptimizationList(basis_type_);

    symmetrizerParser(extractValue<YAML::Node>(custom_node, "symmetrizer"));

    if (extractValue<YAML::Node>(custom_node, "tz_sorter").IsDefined()) {
        optimizations_list_->TzSort();
    }
    if (extractValue<YAML::Node>(custom_node, "tsquared_sorter").IsDefined()) {
        optimizations_list_->TSquaredSort();
    }
    if (extractValue<YAML::Node>(custom_node, "positive_tz_eliminator").IsDefined()) {
        optimizations_list_->EliminatePositiveProjections();
    }
    if (extractValue<YAML::Node>(custom_node, "nonminimal_tz_eliminator").IsDefined()) {
        optimizations_list_->EliminateNonMininalProjections();
    }
    if (extractValue<YAML::Node>(custom_node, "non_abelian_simplifier").IsDefined()) {
        optimizations_list_->NonAbelianSimplify();
    }

    ftlmParser(extractValue<YAML::Node>(custom_node, "ftlm"));

    throw_if_node_is_not_empty(custom_node);
}

void OptimizationsParser::symmetrizerParser(YAML::Node symmetrizer_node) {
    if (!symmetrizer_node.IsDefined()) {
        return;
    }
    // if symmetrizer_node is sequence, it is already correct on its level
    if (!symmetrizer_node.IsSequence()) {
        throw std::invalid_argument("Incorrect format of optimizations::symmetrizer");
    }
    for (auto group_node : symmetrizer_node) {
        groupParser(group_node);
    }
}

void OptimizationsParser::groupParser(YAML::Node group_node) {
    auto group_name = extractValue<group::Group::GroupTypeEnum>(group_node, "group_name");
    std::optional<unsigned int> order = std::nullopt;
    if (group_node["order"].IsDefined()) {
        order = extractValue<unsigned int>(group_node, "order");
    }
    auto generators_vector = extractValue<std::vector<group::Permutation>>(group_node, "generators");

    throw_if_node_is_not_empty(group_node);

    auto group = group::Group({group_name, order}, generators_vector);
    optimizations_list_->Symmetrize(group);
}

void OptimizationsParser::ftlmParser(YAML::Node ftlm_node) {
    if (!ftlm_node.IsDefined()) {
        return;
    }
    auto krylov_subspace_size = extractValue<size_t>(ftlm_node, "krylov_subspace_size");
    auto exact_decomposition_threshold = extractValue<size_t>(ftlm_node, "exact_decomposition_threshold");
    auto number_of_seeds = extractValue<size_t>(ftlm_node, "number_of_seeds");

    throw_if_node_is_not_empty(ftlm_node);

    common::physical_optimization::OptimizationList::FTLMSettings settings;
    settings.exact_decomposition_threshold = exact_decomposition_threshold;
    settings.krylov_subspace_size = krylov_subspace_size;
    settings.number_of_seeds = number_of_seeds;

    optimizations_list_->FTLMApproximate(settings);
}

}  // namespace input