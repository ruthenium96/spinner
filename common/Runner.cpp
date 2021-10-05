#include "Runner.h"
#include "components/NonAbelianSimplifier.h"
#include "components/Symmetrizer.h"
#include "components/TzSorter.h"


runner::Runner::Runner(std::vector<int> mults) : converter_(std::move(mults)), space_(converter_.total_space_size)
{}

void runner::Runner::NonAbelianSimplify() {
    if (history_.number_of_non_simplified_abelian_groups == 0) {
        return;
    }
    if (history_.number_of_non_simplified_abelian_groups != 1) {
        throw std::invalid_argument("Non-Abelian simplification after using of two Non-Abelian Symmetrizers "
                                    "currently is not allowed. Use Non-Abelian simplification twice.");
    }
    NonAbelianSimplifier nonAbelianSimplifier;
    space_ = nonAbelianSimplifier.apply(space_);
    history_.number_of_non_simplified_abelian_groups = 0;
}

void runner::Runner::Symmetrize(Group new_group) {
    // check if user trying to use the same Group for a second time:
    if (std::count(history_.applied_groups.begin(), history_.applied_groups.end(), new_group)) {
        return;
    }

    Symmetrizer symmetrizer(converter_, new_group);
    space_ = symmetrizer.apply(space_);

    if (!IsAbelianGroup(new_group)) {
        ++history_.number_of_non_simplified_abelian_groups;
    }
    history_.applied_groups.emplace_back(std::move(new_group));
}

void runner::Runner::Symmetrize(group::GroupNames group_name, std::vector<Permutation> generators) {
    Group new_group(group_name, std::move(generators));
    Symmetrize(new_group);
}

void runner::Runner::TzSort() {
    // It does not make any sense to use tz_sorter twice.
    if (history_.isTzSorted) {
        return;
    }
    TzSorter tz_sorter(converter_);
    space_ = tz_sorter.apply(space_);
    history_.isTzSorted = true;
}

const Space &runner::Runner::getSpace() const {
    return space_;
}

uint32_t runner::Runner::getTotalSpaceSize() const {
    return converter_.total_space_size;
}

bool runner::Runner::IsAbelianGroup(const Group &group) const noexcept {
    return group.info.group_size == group.info.number_of_representations;
}
