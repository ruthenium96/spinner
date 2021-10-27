//#include "DenseMatrix.h"
//
//#include <armadillo>
//
//struct DenseVector::SubspectrumDataImpl {
//    arma::dvec eigenvalues;
//public:
//    SubspectrumDataImpl() = default;
//};
//
//DenseVector::DenseVector()
//: pImpl{std::make_unique<DenseVector::SubspectrumDataImpl>()} {
//}
//
//std::ostream &operator<<(std::ostream &os, const DenseVector &raw_data) {
//    os << raw_data.pImpl->eigenvalues << std::endl;
//    return os;
//}
//
//DenseVector::~DenseVector() = default;
//DenseVector::DenseVector(DenseVector&&) noexcept = default;
//DenseVector& DenseVector::operator=(DenseVector&&) noexcept = default;
