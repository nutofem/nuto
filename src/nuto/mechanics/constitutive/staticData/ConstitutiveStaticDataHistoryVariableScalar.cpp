#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION


#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataHistoryVariableScalar.h"
#include "nuto/mechanics/elements/ElementEnum.h"



NuTo::ConstitutiveStaticDataHistoryVariableScalar::ConstitutiveStaticDataHistoryVariableScalar():
    ConstitutiveStaticDataBase::ConstitutiveStaticDataBase(), mHistoryVariable(0.0)
{
}

bool NuTo::ConstitutiveStaticDataHistoryVariableScalar::CheckConstitutiveCompatibility(NuTo::Constitutive::eConstitutiveType rConstitutiveType, NuTo::Element::eElementType rElementType) const
{
    if (rConstitutiveType == NuTo::Constitutive::eConstitutiveType::PHASE_FIELD)
    {
        if (rElementType == NuTo::Element::eElementType::CONTINUUMELEMENT)
            return true;
    }
    return false;
}


#ifdef ENABLE_SERIALIZATION
template void NuTo::ConstitutiveStaticDataHistoryVariableScalar::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataHistoryVariableScalar::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataHistoryVariableScalar::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataHistoryVariableScalar::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataHistoryVariableScalar::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstitutiveStaticDataHistoryVariableScalar::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstitutiveStaticDataHistoryVariableScalar::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstitutiveStaticDataHistoryVariableScalar" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveStaticDataBase)
       & BOOST_SERIALIZATION_NVP(mHistoryVariable);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstitutiveStaticDataHistoryVariableScalar" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstitutiveStaticDataHistoryVariableScalar)
#endif // ENABLE_SERIALIZATION
