#pragma once

#include "components/symmetrizer.h"
#include "components/tz_sorter.h"
#include "entities/Space.h"

namespace runner {
class Runner {
  public:
    Runner(std::vector<int> mults = {3, 3, 3})
        : converter_(std::move(mults)), 
          space_(converter_.total_space_size), 
          tz_sorter_(converter_),
          group_first_(group::S3, {{1, 2, 0}, {0, 2, 1}}),
          symmetrizer_first(converter_, group_first_) {
    }

    void TzSort() {
        // It does not make any sense to use tz_sorter twice.
        if (history_.isTzSorted) {
            return;
        }
        space_ = tz_sorter_.apply(space_);
        history_.isTzSorted = true;
    }

    void Symmetrize() {
        space_ = symmetrizer_first.apply(space_);
        history_.isC2Symmetrized = true;
    }
    
    const Space& getSpace() const {
        return space_;
    }

    void PrintSpace() const {
        std::cout << space_;
    }

  private:
    struct History {
        bool isTzSorted = false;
        bool isC2Symmetrized = false;
    };
    spaces::LexicographicIndexConverter converter_;
    Space space_;
    Tz_Sorter tz_sorter_;
    Group group_first_;
    Symmetrizer symmetrizer_first;
    History history_;
};
} // namespace runner