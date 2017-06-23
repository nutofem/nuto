#pragma once

#include <vector>

namespace NuTo {
namespace Math
{

class Polynomial
{
public:
    //! value of the k-th derivative of legendre polynomial
    //! of the first kind of order n
    static double legendre_p_deriv(int n, double x , int k);

    //! roots of legendre polynomial
    //! of the first kind of order n
    static std::vector<double> legendre_roots(int n);

    //! roots of derivative of legendre polynomial
    //! of the first kind of order n
    static std::vector<double> legendre_deriv_roots(int n);
};


} // namespace Math
} // namespace NuTo


