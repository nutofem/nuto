// $Id: ConstitutiveStaticDataBase.cpp 102 2009-11-11 10:47:23Z eckardt4 $


#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataMisesPlasticity.h"
#include "nuto/mechanics/MechanicsException.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION


#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"

//! @brief constructor
template<int TDim>
NuTo::ConstitutiveStaticDataMisesPlasticity<TDim>::ConstitutiveStaticDataMisesPlasticity() :
        ConstitutiveStaticDataBase(), mEpsilonPEq(0)
{
    mEpsilonP.setZero();
    mSigmaB.setZero();
}

//! @brief check, if the static data is compatible with a given element and a given constitutive model
template<int TDim>
bool NuTo::ConstitutiveStaticDataMisesPlasticity<TDim>::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const
{
    if (rConstitutiveType == NuTo::Constitutive::eConstitutiveType::MISES_PLASTICITY_ENGINEERING_STRESS)
    {
        if (rElementType == Element::eElementType::CONTINUUMELEMENT)
            return true;

    }

    return false;
}

namespace NuTo
{
template<>
ConstitutiveStaticDataMisesPlasticity<1>* ConstitutiveStaticDataMisesPlasticity<1>::AsConstitutiveStaticDataMisesPlasticity1D()
{
    return this;
}
template<>
const ConstitutiveStaticDataMisesPlasticity<1>* ConstitutiveStaticDataMisesPlasticity<1>::AsConstitutiveStaticDataMisesPlasticity1D() const
{
    return this;
}
template<>
ConstitutiveStaticDataMisesPlasticity<2>* ConstitutiveStaticDataMisesPlasticity<2>::AsConstitutiveStaticDataMisesPlasticity2D()
{
    return this;
}
template<>
const ConstitutiveStaticDataMisesPlasticity<2>* ConstitutiveStaticDataMisesPlasticity<2>::AsConstitutiveStaticDataMisesPlasticity2D() const
{
    return this;
}
template<>
ConstitutiveStaticDataMisesPlasticity<3>* ConstitutiveStaticDataMisesPlasticity<3>::AsConstitutiveStaticDataMisesPlasticity3D()
{
    return this;
}
template<>
const ConstitutiveStaticDataMisesPlasticity<3>* ConstitutiveStaticDataMisesPlasticity<3>::AsConstitutiveStaticDataMisesPlasticity3D() const
{
    return this;
}



}
  // namespace NuTo

template class NuTo::ConstitutiveStaticDataMisesPlasticity<1>;
template class NuTo::ConstitutiveStaticDataMisesPlasticity<2>;
template class NuTo::ConstitutiveStaticDataMisesPlasticity<3>;


#ifdef ENABLE_SERIALIZATION

template void NuTo::ConstitutiveStaticDataMisesPlasticity<1>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<2>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<3>::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<1>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<2>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<3>::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<1>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<2>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<3>::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<1>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<2>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<3>::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<1>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<2>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<3>::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<1>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<2>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataMisesPlasticity<3>::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<int TDim>
template<class Archive> void NuTo::ConstitutiveStaticDataMisesPlasticity<TDim>::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataMisesPlasticity" << std::endl;
    #endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
    & BOOST_SERIALIZATION_NVP(mEpsilonPEq)
    & BOOST_SERIALIZATION_NVP(mEpsilonP)
    & BOOST_SERIALIZATION_NVP(mSigmaB);
    #ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataMisesPlasticity" << std::endl;
    #endif
}

BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ConstitutiveStaticDataMisesPlasticity<1>)), "ConstitutiveStaticDataMisesPlasticity_1")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ConstitutiveStaticDataMisesPlasticity<2>)), "ConstitutiveStaticDataMisesPlasticity_2")
BOOST_CLASS_EXPORT_GUID(BOOST_IDENTITY_TYPE((NuTo::ConstitutiveStaticDataMisesPlasticity<3>)), "ConstitutiveStaticDataMisesPlasticity_3")
#endif // ENABLE_SERIALIZATION
