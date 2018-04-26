#pragma once

#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_resize.hpp>

namespace Eigen
{

Eigen::VectorXd operator+(const Eigen::VectorXd& v, double scalar)
{
    return v + Eigen::VectorXd::Constant(v.rows(), scalar);
}
Eigen::VectorXd operator+(double scalar, const Eigen::VectorXd& v)
{
    return v + Eigen::VectorXd::Constant(v.rows(), scalar);
}

template <typename D1, typename D2>
inline const typename Eigen::CwiseBinaryOp<
        typename Eigen::internal::scalar_quotient_op<typename Eigen::internal::traits<D1>::Scalar>, const D1, const D2>
operator/(const Eigen::MatrixBase<D1>& x1, const Eigen::MatrixBase<D2>& x2)
{
    return x1.cwiseQuotient(x2);
}

template <typename D>
inline const typename Eigen::CwiseUnaryOp<
        typename Eigen::internal::scalar_abs_op<typename Eigen::internal::traits<D>::Scalar>, const D>
abs(const Eigen::MatrixBase<D>& m)
{
    return m.cwiseAbs();
}

} // end Eigen namespace

namespace boost
{
namespace numeric
{
namespace odeint
{

// needed for steppers with error control
template <class Derivative>
struct vector_space_norm_inf<Eigen::MatrixBase<Derivative>>
{
    typedef double result_type;
    double operator()(const Eigen::MatrixBase<Derivative>& m) const
    {
        return m.template lpNorm<Eigen::Infinity>();
    }
};

}
}
} // end boost::numeric::odeint namespace
