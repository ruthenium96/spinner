#include "UnitarySparseMatrix.h"
#include <armadillo>

struct NewBasisDecomposition::Impl {
    arma::sp_mat sparse_basis;
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
    return std::make_unique<IteratorImpl>(pImpl->sparse_basis.begin_col(index_of_vector),
                                          pImpl->sparse_basis.end_col(index_of_vector));
}

std::ostream &operator<<(std::ostream &os, const NewBasisDecomposition &decomposition) {
    os << decomposition.pImpl->sparse_basis << std::endl;
    return os;
}

uint32_t NewBasisDecomposition::size() const {
    return pImpl->sparse_basis.n_cols;
}

bool NewBasisDecomposition::empty() const {
    return pImpl->sparse_basis.n_cols == 0;
}

bool NewBasisDecomposition::vempty(uint32_t index_of_vector) const {
    return pImpl->sparse_basis.col(index_of_vector).n_nonzero == 0;
}

void NewBasisDecomposition::clear() {
    pImpl->sparse_basis.reset();
}

void NewBasisDecomposition::move_vector_from(uint32_t i, NewBasisDecomposition& decomposition_from) {
    tensor_size = decomposition_from.tensor_size;
    pImpl->sparse_basis.resize(tensor_size, pImpl->sparse_basis.n_cols + 1);
    pImpl->sparse_basis.col(pImpl->sparse_basis.n_cols - 1) = decomposition_from.pImpl->sparse_basis.col(i);
    // TODO: Is it good idea?
//    decomposition_from.pImpl->sparse_basis.col(i).zeros();
}

void NewBasisDecomposition::move_all_from(NewBasisDecomposition& decomposition_from) {
    tensor_size = decomposition_from.tensor_size;
    pImpl->sparse_basis.resize(tensor_size, pImpl->sparse_basis.n_cols + decomposition_from.pImpl->sparse_basis.n_cols);
    pImpl->sparse_basis.cols(pImpl->sparse_basis.n_cols - decomposition_from.pImpl->sparse_basis.n_cols, pImpl->sparse_basis.n_cols - 1) = decomposition_from.pImpl->sparse_basis;
    decomposition_from.clear();
}

void NewBasisDecomposition::copy_vector_from(uint32_t i, const NewBasisDecomposition& decomposition_from) {
    tensor_size = decomposition_from.tensor_size;
    pImpl->sparse_basis.resize(tensor_size, pImpl->sparse_basis.n_cols + 1);
    pImpl->sparse_basis.col(pImpl->sparse_basis.n_cols - 1) = decomposition_from.pImpl->sparse_basis.col(i);
}

void NewBasisDecomposition::copy_all_from(const NewBasisDecomposition& decomposition_from) {
    tensor_size = decomposition_from.tensor_size;
    pImpl->sparse_basis.resize(tensor_size, pImpl->sparse_basis.n_cols + decomposition_from.pImpl->sparse_basis.n_cols);
    pImpl->sparse_basis.cols(pImpl->sparse_basis.n_cols - decomposition_from.pImpl->sparse_basis.n_cols, pImpl->sparse_basis.n_cols - 1) = decomposition_from.pImpl->sparse_basis;
}

void NewBasisDecomposition::add_to_position(double value, uint32_t i, uint32_t j) {
    pImpl->sparse_basis(j, i) += value;
}

double NewBasisDecomposition::operator()(uint32_t i, uint32_t j) const {
    return pImpl->sparse_basis(j, i);
}

void NewBasisDecomposition::resize(uint32_t new_size) {
    pImpl->sparse_basis.resize(tensor_size, new_size);
}

bool NewBasisDecomposition::is_zero(uint32_t i, uint32_t j) const {
    // TODO: epsilon
    return std::abs(pImpl->sparse_basis(j, i)) < 0.001;
}

void NewBasisDecomposition::erase_if_zero() {
    for (auto && el : pImpl->sparse_basis) {
        // TODO: epsilon
        if (std::abs(el) < 0.001) {
            el = 0.0;
        }
    }
}