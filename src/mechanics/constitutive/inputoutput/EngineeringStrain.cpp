#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"

template <int TDim>
NuTo::EngineeringStrain<TDim>::EngineeringStrain(std::initializer_list<double> rList)
{
    assert(rList.size() == ConstitutiveIOBase::GetVoigtDim(TDim));
    for (auto iterator = rList.begin(); iterator != rList.end(); ++iterator)
    {
        int position = std::distance(rList.begin(), iterator);
        (*this)[position] = *iterator;
    }
}


template<int TDim>
std::unique_ptr<NuTo::ConstitutiveIOBase> NuTo::EngineeringStrain<TDim>::clone()
{
    return std::make_unique<EngineeringStrain<TDim>>(*this);
}

template<>
template<int U, typename std::enable_if<U == 3, int>::type>
double NuTo::EngineeringStrain<3>::InvariantI1() const
{
    return (*this).segment<3>(0).sum(); // sum of the first 3 entries
}

template<>
template<int U, typename std::enable_if<U == 3, int>::type>
double NuTo::EngineeringStrain<3>::InvariantI2() const
{
    const Eigen::Matrix<double, 6, 1>& e = *this;
    // e_xx = v[0]
    // e_yy = v[1]
    // e_zz = v[2]
    // e_yz = v[3] / 2
    // e_zx = v[4] / 2
    // e_xy = v[5] / 2
    return e[0] * e[1] + e[1] * e[2] + e[2] * e[0] - 0.25 * (e[3] * e[3] + e[4] * e[4] + e[5] * e[5]);
}

template<>
template<int U, typename std::enable_if<U == 3, int>::type>
double NuTo::EngineeringStrain<3>::InvariantI3() const
{
    const Eigen::Matrix<double, 6, 1>& v = *this;
    // e_xx = v[0]
    // e_yy = v[1]
    // e_zz = v[2]
    // e_yz = v[3] / 2
    // e_zx = v[4] / 2
    // e_xy = v[5] / 2
    return v[0] * v[1] * v[2] + 0.25 * v[3] * v[4] * v[5] - .25 * (v[5] * v[5] * v[2]  +  v[3] * v[3] * v[0]  +  v[4] * v[4] * v[1]);
}

template<>
template<int U, typename std::enable_if<U == 3, int>::type>
double NuTo::EngineeringStrain<3>::InvariantJ2() const
{
    const Eigen::Matrix<double, 6, 1>& v = *this;
    const double eps_xx = v[0];
    const double eps_yy = v[1];
    const double eps_zz = v[2];
    const double eps_yz = .5*v[3];
    const double eps_zx = .5*v[4];
    const double eps_xy = .5*v[5];


    return 1./6. * ( (eps_xx - eps_yy)*(eps_xx - eps_yy) + (eps_yy - eps_zz)*(eps_yy - eps_zz) + (eps_zz - eps_xx)*(eps_zz - eps_xx) )
            + eps_xy*eps_xy +eps_yz*eps_yz + eps_zx*eps_zx;
}


template<>
template<int U, typename std::enable_if<U == 3, int>::type>
NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<3>::Deviatoric() const
{
    EngineeringStrain<3> v = *this;
    double I1_3 = v.InvariantI1()/3.;
    v[0] -= I1_3;
    v[1] -= I1_3;
    v[2] -= I1_3;
    return v;
}

namespace NuTo
{



template<>
NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<1>::As3D(double rNu, ePlaneState) const
{
    EngineeringStrain<3> strain3D;
    double strain = (*this)[0];
    strain3D[0] = strain;
    strain3D[1] = -rNu * strain;
    strain3D[2] = -rNu * strain;
    strain3D[3] = 0.;
    strain3D[4] = 0.;
    strain3D[5] = 0;
    return strain3D;
}

template<>
NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<2>::As3D(double rNu, ePlaneState planeState) const
{
    EngineeringStrain<3> strain3D;
    const EngineeringStrain<2>& v = *this;
    switch (planeState)
    {
    case ePlaneState::PLANE_STRAIN:
        strain3D[0] = v[0];
        strain3D[1] = v[1];
        strain3D[2] = 0;
        strain3D[3] = 0.;
        strain3D[4] = 0.;
        strain3D[5] = v[2];
        break;
    case ePlaneState::PLANE_STRESS:
        strain3D[0] = v[0];
        strain3D[1] = v[1];
        strain3D[2] = rNu / (rNu - 1.) * (v[0] + v[1]);
        strain3D[3] = 0.;
        strain3D[4] = 0.;
        strain3D[5] = v[2];
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, 
                "Specify section behavior, either PLANE_STRAIN or PLANE_STRESS");
    }
    return strain3D;
}

template<>
NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<3>::As3D(double, ePlaneState) const
{
    return *this;
}
}  // namespace NuTo



template double NuTo::EngineeringStrain<3>::InvariantI1() const;
template double NuTo::EngineeringStrain<3>::InvariantI2() const;
template double NuTo::EngineeringStrain<3>::InvariantI3() const;
template double NuTo::EngineeringStrain<3>::InvariantJ2() const;

template NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<3>::Deviatoric() const;


template NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<1>::As3D(double, ePlaneState) const;
template NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<2>::As3D(double, ePlaneState) const;
template NuTo::EngineeringStrain<3> NuTo::EngineeringStrain<3>::As3D(double, ePlaneState) const;


template class NuTo::EngineeringStrain<1>;
template class NuTo::EngineeringStrain<2>;
template class NuTo::EngineeringStrain<3>;
