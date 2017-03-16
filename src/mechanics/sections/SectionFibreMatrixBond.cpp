#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include <iostream>
#include "mechanics/sections/SectionFibreMatrixBond.h"

using namespace NuTo;


SectionFibreMatrixBond::SectionFibreMatrixBond(double circumference):
mCircumference(circumference)
{
}


std::shared_ptr<SectionFibreMatrixBond> SectionFibreMatrixBond::Create(double circumference)
{
    return std::shared_ptr<SectionFibreMatrixBond>(new SectionFibreMatrixBond(circumference));
}


double NuTo::SectionFibreMatrixBond::GetCircumference() const
{
    return mCircumference;
}


void NuTo::SectionFibreMatrixBond::Info() const
{
    std::cout << "    Fibre matric bond section with circumference: " << mCircumference << std::endl;
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
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Section)
       & BOOST_SERIALIZATION_NVP(mCircumference);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SectionFibreMatrixBond" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionFibreMatrixBond)
#endif // ENABLE_SERIALIZATION
