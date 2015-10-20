#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataBondStressSlip.h"
#include "nuto/mechanics/MechanicsException.h"

NuTo::ConstitutiveStaticDataBondStressSlip::ConstitutiveStaticDataBondStressSlip() :
        NuTo::ConstitutiveStaticDataBase::ConstitutiveStaticDataBase(), mSlip(0)
{
}

NuTo::ConstitutiveStaticDataBondStressSlip& NuTo::ConstitutiveStaticDataBondStressSlip::operator=(NuTo::ConstitutiveStaticDataBondStressSlip const& rOther)
{
    NuTo::ConstitutiveStaticDataBase::operator=(rOther);
    mSlip = rOther.mSlip;
    return *this;
}

void NuTo::ConstitutiveStaticDataBondStressSlip::SetSlip(double rSlip)
{
    mSlip = rSlip;
}

const double NuTo::ConstitutiveStaticDataBondStressSlip::GetSlip() const
{
    return mSlip;
}

NuTo::ConstitutiveStaticDataBondStressSlip* NuTo::ConstitutiveStaticDataBondStressSlip::AsBondStressSlip()
{
    return this;
}

const NuTo::ConstitutiveStaticDataBondStressSlip* NuTo::ConstitutiveStaticDataBondStressSlip::AsBondStressSlip() const
{
    return this;
}

NuTo::ConstitutiveStaticDataBondStressSlip* NuTo::ConstitutiveStaticDataBondStressSlip::Clone()const
{
    return new ConstitutiveStaticDataBondStressSlip(*this);
}

bool NuTo::ConstitutiveStaticDataBondStressSlip::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const
{
    if (rConstitutiveType == NuTo::Constitutive::INTERFACE_GOODMAN)
    {
        if (rElementType == NuTo::Element::ELEMENT2DINTERFACE)
            return true;
        else
            return false;
    } else
        return false;
}

