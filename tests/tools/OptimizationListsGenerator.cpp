#include "OptimizationListsGenerator.h"

std::vector<common::physical_optimization::OptimizationList> generate_all_optimization_lists(const std::vector<group::Group> groups) {
    std::vector<common::physical_optimization::OptimizationList> answer = 
        {
            common::physical_optimization::OptimizationList(common::physical_optimization::OptimizationList::LEX),
            common::physical_optimization::OptimizationList(common::physical_optimization::OptimizationList::ITO)
        };
    size_t last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        try {
            common::physical_optimization::OptimizationList copy = answer[i];
            copy.TzSort();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }
    last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        try {
            common::physical_optimization::OptimizationList copy = answer[i];
            copy.EliminatePositiveProjections();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }
    for (const auto& group : groups) {
        last_size = answer.size();
        for (int i = 0; i < last_size; ++i) {
            try {
                common::physical_optimization::OptimizationList copy = answer[i];
                copy.Symmetrize(group);
                answer.push_back(copy);
            } catch (std::invalid_argument) {}
        }
    }

    last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        try {
            common::physical_optimization::OptimizationList copy = answer[i];
            copy.NonAbelianSimplify();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }

    last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        if (answer[i].isLexBasis()) {
            continue;
        }
        try {
            common::physical_optimization::OptimizationList copy = answer[i];
            copy.TSquaredSort();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }
    last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        if (answer[i].isLexBasis()) {
            continue;
        }
        try {
            common::physical_optimization::OptimizationList copy = answer[i];
            copy.EliminateNonMininalProjections();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }

    return answer;
}

