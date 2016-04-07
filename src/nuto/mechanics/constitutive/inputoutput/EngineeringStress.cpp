#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"



namespace NuTo
{

template<>
double EngineeringStress<1>::GetVonMisesStress() const
{
    return (*this)[0];
}

template<>
double EngineeringStress<2>::GetVonMisesStress() const
{
    throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] implement me!");
}

template<>
double EngineeringStress<3>::GetVonMisesStress() const
{
    throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] implement me!");
}


} /* namespace NuTo */


template class NuTo::EngineeringStress<1>;
template class NuTo::EngineeringStress<2>;
template class NuTo::EngineeringStress<3>;
