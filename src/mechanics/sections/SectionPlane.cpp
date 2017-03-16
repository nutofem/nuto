#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/MechanicsException.h"
#include "mechanics/sections/SectionPlane.h"

using namespace NuTo;

SectionPlane::SectionPlane(double thickness, bool isPlaneStrain)
    : mThickness(thickness)
    , mIsPlaneStrain(isPlaneStrain)
{
}


std::shared_ptr<SectionPlane> SectionPlane::Create(double thickness, bool isPlaneStrain)
{
    return std::shared_ptr<SectionPlane>(new SectionPlane(thickness, isPlaneStrain));
}


double SectionPlane::GetThickness() const
{
    return mThickness;
}


void SectionPlane::Info() const
{
    std::cout << "    Plane section with thickness: " << mThickness << std::endl;
    if (mIsPlaneStrain)
        std::cout << "    Section type is plane strain." << std::endl;
    else
        std::cout << "    Section type is plane stress." << std::endl;
}


bool SectionPlane::IsPlaneStrain() const
{
    return mIsPlaneStrain;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void SectionPlane::serialize(boost::archive::binary_oarchive& ar, const unsigned int version);
template void SectionPlane::serialize(boost::archive::xml_oarchive& ar, const unsigned int version);
template void SectionPlane::serialize(boost::archive::text_oarchive& ar, const unsigned int version);
template void SectionPlane::serialize(boost::archive::binary_iarchive& ar, const unsigned int version);
template void SectionPlane::serialize(boost::archive::xml_iarchive& ar, const unsigned int version);
template void SectionPlane::serialize(boost::archive::text_iarchive& ar, const unsigned int version);
template <class Archive>
void SectionPlane::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionPlane" << std::endl;
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(Section) & BOOST_SERIALIZATION_NVP(mThickness) &
            BOOST_SERIALIZATION_NVP(mSectionType);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SectionPlane" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(SectionPlane)
#endif // ENABLE_SERIALIZATION
