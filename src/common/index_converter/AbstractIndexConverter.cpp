#include "AbstractIndexConverter.h"
#include <numeric>

index_converter::AbstractIndexConverter::AbstractIndexConverter(
	std::vector<spin_algebra::Multiplicity> mults) : mults_(mults) {

    // It is the size of our projection-based space.
	total_space_size_ = std::accumulate(mults_.begin(), 
										mults_.end(), 
										1, 
										std::multiplies<spin_algebra::Multiplicity>());

	spins_.resize(mults_.size());
    for (size_t i = 0; i < mults.size(); ++i) {
        spins_[i] = (mults_[i] - 1) / 2.0;
    }

    // We want to get 2T + 1 (projections are counted from zero to multiplicity),
    // where T = sum_{1}^{N} S_i. So 2T + 1 = sum_{1}^{N} (2S_i + 1) - N + 1.
    max_ntz_proj_= std::accumulate(get_mults().begin(), get_mults().end(), 1 - get_mults().size());

}

const std::vector<spin_algebra::Multiplicity>&
index_converter::AbstractIndexConverter::get_mults() const {
    return mults_;
}

const std::vector<double>& index_converter::AbstractIndexConverter::get_spins() const {
	return spins_;
}

uint32_t index_converter::AbstractIndexConverter::get_total_space_size() const {
    return total_space_size_;
}

uint32_t index_converter::AbstractIndexConverter::get_max_ntz_proj() const {
    return max_ntz_proj_;
}
