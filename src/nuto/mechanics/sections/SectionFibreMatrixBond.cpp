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
#include "nuto/mechanics/sections/SectionFibreMatrixBond.h"

// constructor
NuTo::SectionFibreMatrixBond::SectionFibreMatrixBond():
mCircumference(0.0)
{
}

// set section circumference
void NuTo::SectionFibreMatrixBond::SetCircumference(double rCircumference)
{
    if (rCircumference <= 0.0)
    {
        throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + ":\t circumference must be greater than 0.") ;
    }
    mCircumference = rCircumference;
}

// get section circumference
double NuTo::SectionFibreMatrixBond::GetCircumference() const
{
    return mCircumference;
}

// get section type
NuTo::Section::eSectionType NuTo::SectionFibreMatrixBond::GetType() const
{
    return NuTo::Section::FIBRE_MATRIX_BOND;
}

// info routine
void NuTo::SectionFibreMatrixBond::Info(unsigned short rVerboseLevel) const
{
    this->SectionBase::Info(rVerboseLevel);
    std::cout << "    section circumference: " << mCircumference << std::endl;

}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::SectionFibreMatrixBond::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SectionFibreMatrixBond::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SectionFibreMatrixBond::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::SectionFibreMatrixBond::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SectionFibreMatrixBond::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SectionFibreMatrixBond::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::SectionFibreMatrixBond::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionFibreMatrixBond" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(SectionBase)
       & BOOST_SERIALIZATION_NVP(mCircumference);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SectionFibreMatrixBond" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionFibreMatrixBond)
#endif // ENABLE_SERIALIZATION
