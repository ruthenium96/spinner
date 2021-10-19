#include "NewBasisDecomposition.h"
#include <armadillo>

struct NewBasisDecomposition::Impl {
    std::vector<arma::sp_vec> sparse_vectors;
public:
    Impl() = default;
};

NewBasisDecomposition::NewBasisDecomposition() : pImpl{std::make_unique<Impl>()} {}
NewBasisDecomposition::~NewBasisDecomposition() = default;
NewBasisDecomposition::NewBasisDecomposition(NewBasisDecomposition&&) noexcept = default;
NewBasisDecomposition& NewBasisDecomposition::operator=(NewBasisDecomposition&&) noexcept = default;


struct IteratorImpl : public NewBasisDecomposition::Iterator {

    arma::SpMat<double>::col_iterator iter;
    const arma::SpMat<double>::col_iterator end;

    IteratorImpl(arma::SpMat<double>::col_iterator iter1,
                 arma::SpMat<double>::col_iterator iter2): iter(iter1), end(iter2){
    }

    [[nodiscard]] bool hasNext() const override {
        return iter != end;
    }

    IndexValueItem getNext() override {
        auto pair = iter;
        ++iter;
        return {static_cast<uint32_t>(pair.row()), *pair};
    }

    ~IteratorImpl() = default;
};


std::unique_ptr<NewBasisDecomposition::Iterator> NewBasisDecomposition::GetNewIterator(size_t index_of_vector) const {
    return std::make_unique<IteratorImpl>(pImpl->sparse_vectors[index_of_vector].begin_col(0),
                                          pImpl->sparse_vectors[index_of_vector].end_col(0));
}

std::ostream &operator<<(std::ostream &os, const NewBasisDecomposition &decomposition) {
    for (const auto& v : decomposition.pImpl->sparse_vectors) {
        os << v << std::endl;
    }
    return os;
}

uint32_t NewBasisDecomposition::size() const {
    return pImpl->sparse_vectors.size();
}

bool NewBasisDecomposition::empty() const {
    return pImpl->sparse_vectors.empty();
}

bool NewBasisDecomposition::vempty(uint32_t index_of_vector) const {
    return pImpl->sparse_vectors[index_of_vector].n_nonzero == 0;
}

void NewBasisDecomposition::clear() {
    pImpl->sparse_vectors.clear();
}


void NewBasisDecomposition::move_vector_from(uint32_t i, NewBasisDecomposition& decomposition_from) {
    pImpl->sparse_vectors.emplace_back(std::move(decomposition_from.pImpl->sparse_vectors[i]));
}

void NewBasisDecomposition::move_all_from(NewBasisDecomposition& decomposition_from) {
    for (uint32_t i = 0; i < decomposition_from.size(); ++i) {
        move_vector_from(i, decomposition_from);
    }
}

void NewBasisDecomposition::copy_vector_from(uint32_t i, const NewBasisDecomposition& decomposition_from) {
    pImpl->sparse_vectors.emplace_back(decomposition_from.pImpl->sparse_vectors[i]);
}

void NewBasisDecomposition::copy_all_from(const NewBasisDecomposition& decomposition_from) {
    for (uint32_t i = 0; i < decomposition_from.size(); ++i) {
        copy_vector_from(i, decomposition_from);
    }
}

void NewBasisDecomposition::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->sparse_vectors[i][j] += value;
}

double NewBasisDecomposition::operator()(uint32_t i, uint32_t j) const {
    return pImpl->sparse_vectors[i][j];
}

void NewBasisDecomposition::resize(uint32_t new_size) {
    pImpl->sparse_vectors.resize(new_size, arma::SpCol<double>(tensor_size));
}

bool NewBasisDecomposition::is_zero(uint32_t i, uint32_t j) const {
    // TODO: epsilon
    return std::abs(pImpl->sparse_vectors[i][j]) < 0.001;
}

void NewBasisDecomposition::erase_if_zero() {
    for (auto & sp_v : pImpl->sparse_vectors) {
        // TODO: epsilon
        for (auto i = sp_v.begin(); i != sp_v.end(); ++i) {
            if (std::abs(*i) < 0.001) {
                *i = 0;
            }
        }
    }
}

// TODO: test if it is really work
void NewBasisDecomposition::normalize() {
    for (auto& v : pImpl->sparse_vectors) {
        double norm = arma::norm(v);
        for (auto&& el : v) {
            el /= norm;
        }
    }
}
