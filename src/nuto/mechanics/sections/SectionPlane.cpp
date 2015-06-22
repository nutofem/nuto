// $Id$
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/sections/SectionPlane.h"

// constructor
NuTo::SectionPlane::SectionPlane(Section::eSectionType rSectionType)
{
    this->mThickness = 0.0;
    if (rSectionType==Section::PLANE_STRAIN || rSectionType==Section::PLANE_STRESS)
    	this->mSectionType = rSectionType;
    else
    {
    	throw NuTo::MechanicsException("[NuTo::SectionPlane::SectionPlane] section for plane elements is either plane strain or plane stress.") ;
    }
}

// set section thickness
void NuTo::SectionPlane::SetThickness(double rThickness)
{
    if (rThickness <= 0.0)
    {
        throw NuTo::MechanicsException("[NuTo::SectionPlane::SetThickness] section thickness must be greater than 0.") ;
    }
    this->mThickness = rThickness;
}

// get section thickness
double NuTo::SectionPlane::GetThickness() const
{
    return this->mThickness;
}

// get section type
NuTo::Section::eSectionType NuTo::SectionPlane::GetType() const
{
    return mSectionType;
}

// info routine
void NuTo::SectionPlane::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
    std::cout << "    section type: 2D" << std::endl;
    std::cout << "    section thickness: " << this->mThickness << std::endl;
    if (mSectionType==Section::PLANE_STRAIN)
        std::cout << "    section stress state: plane strain" << std::endl;
    else if (mSectionType==Section::PLANE_STRESS)
        std::cout << "    section stress state: plane stress" << std::endl;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::SectionPlane::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SectionPlane::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SectionPlane::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::SectionPlane::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SectionPlane::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SectionPlane::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::SectionPlane::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionPlane" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SectionBase)
       & BOOST_SERIALIZATION_NVP(mThickness)
       & BOOST_SERIALIZATION_NVP(mSectionType);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SectionPlane" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionPlane)
#endif // ENABLE_SERIALIZATION
