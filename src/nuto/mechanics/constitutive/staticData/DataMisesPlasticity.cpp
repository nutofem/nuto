#include "nuto/mechanics/constitutive/staticData/DataMisesPlasticity.h"

using namespace NuTo::Constitutive::StaticData;

template<int TDim>
double DataMisesPlasticity<TDim>::GetEquivalentPlasticStrain() const
{
    return mEpsilonPEq;
}


template<int TDim>
NuTo::EngineeringStrain<TDim> DataMisesPlasticity<TDim>::GetPlasticStrain() const
{
    return mEpsilonP;
}


template<int TDim>
NuTo::EngineeringStress<TDim> DataMisesPlasticity<TDim>::GetBackStress() const
{
    return mSigmaB;
}


template<int TDim>
void DataMisesPlasticity<TDim>::SetEquivalentPlasticStrain(double newEpsilonPEq)
{
    mEpsilonPEq = newEpsilonPEq;
}


template<int TDim>
void DataMisesPlasticity<TDim>::SetPlasticStrain(EngineeringStrain<TDim> newEpsilonP)
{
    mEpsilonP = newEpsilonP;
}

  
template<int TDim>
void DataMisesPlasticity<TDim>::SetBackStress(EngineeringStress<TDim> newBackStress)
{
    mSigmaB = newBackStress;
}

namespace NuTo
{
namespace Constitutive
{
namespace StaticData
{
    template class DataMisesPlasticity<1>;
    template class DataMisesPlasticity<2>;
    template class DataMisesPlasticity<3>;
}
}
}
