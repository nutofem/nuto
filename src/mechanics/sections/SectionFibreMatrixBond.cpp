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


double SectionFibreMatrixBond::GetCircumference() const
{
    return mCircumference;
}


void SectionFibreMatrixBond::Info(std::ostream& out) const
{
    out << "    Fibre matrix bond section with circumference: " << mCircumference << "\n";
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void SectionFibreMatrixBond::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void SectionFibreMatrixBond::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void SectionFibreMatrixBond::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void SectionFibreMatrixBond::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void SectionFibreMatrixBond::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void SectionFibreMatrixBond::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void SectionFibreMatrixBond::serialize(Archive & ar, const unsigned int version)
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
BOOST_CLASS_EXPORT_IMPLEMENT(SectionFibreMatrixBond)
#endif // ENABLE_SERIALIZATION
