#include "nuto/mechanics/constitutive/staticData/DataMisesPlasticity.h"

using namespace NuTo::Constitutive::StaticData;

double DataMisesPlasticity::GetEquivalentPlasticStrain() const
{
    return mEpsilonPEq;
}


EngineeringStrain<TDim> DataMisesPlasticity::GetPlasticStrain() const
{
    return mEpsilonP;
}


EngineeringStress<TDim> DataMisesPlasticity::GetBackStress() const
{
    return mSigmaB;
}


void SetEquivalentPlasticStrain(double newEpsilonPEq)
{
    mEpsilonPEq = newEpsilonPEq;
}


void SetPlasticStrain(EngineeringStrain<TDim> newEpsilonP)
{
    mEpsilonP = newEpsilonP;
}

  
void SetBackStress(EngineeringStress<TDim> newBackStress)
{
    mSigmaB = newBackStress;
}
