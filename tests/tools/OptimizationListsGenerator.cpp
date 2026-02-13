#include "OptimizationListsGenerator.h"

#include <stdexcept>

using namespace spinner::common::physical_optimization;

std::vector<OptimizationList> generate_all_optimization_lists(const std::vector<spinner::group::Group> groups) {
    std::vector<OptimizationList> answer = 
        {
            OptimizationList(OptimizationList::LEX),
            OptimizationList(OptimizationList::ITO)
        };
    size_t last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        try {
            OptimizationList copy = answer[i];
            copy.TzSort();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }
    last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        try {
            OptimizationList copy = answer[i];
            copy.EliminatePositiveProjections();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }
    for (const auto& group : groups) {
        last_size = answer.size();
        for (int i = 0; i < last_size; ++i) {
            try {
                OptimizationList copy = answer[i];
                copy.Symmetrize(group);
                answer.push_back(copy);
            } catch (std::invalid_argument) {}
        }
    }

    last_size = answer.size();
    for (int i = 0; i < last_size; ++i) {
        try {
            OptimizationList copy = answer[i];
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
            OptimizationList copy = answer[i];
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
            OptimizationList copy = answer[i];
            copy.EliminateNonMininalProjections();
            answer.push_back(copy);
        } catch (std::invalid_argument) {}
    }

    return answer;
}

