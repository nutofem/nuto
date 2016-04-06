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
    if (rConstitutiveType == NuTo::Constitutive::FIBRE_MATRIX_BOND_STRESS_SLIP)
    {
        if (rElementType == NuTo::Element::ELEMENT2DINTERFACE)
            return true;
        else
            return false;
    } else
        return false;
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template void NuTo::ConstitutiveStaticDataBondStressSlip::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBondStressSlip::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBondStressSlip::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBondStressSlip::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBondStressSlip::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataBondStressSlip::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataBondStressSlip::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataBondStressSlip" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mSlip);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataBondStressSlip" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataBondStressSlip)
#endif // ENABLE_SERIALIZATION
