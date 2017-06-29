#include <iostream>
#include "base/Exception.h"
#include "mechanics/sections/SectionTruss.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/vector.hpp>
#endif

using namespace NuTo;

SectionTruss::SectionTruss(double area)
    : mArea(area)
{
}

std::shared_ptr<SectionTruss> SectionTruss::Create(double area)
{
    return std::shared_ptr<SectionTruss>(new SectionTruss(area));
}


double SectionTruss::GetArea(double) const
{
    return mArea;
}


void SectionTruss::Info(std::ostream& out) const
{
    out << "    Truss section with area " << mArea << "\n";
}

#ifdef ENABLE_SERIALIZATION
template <class Archive>
void NuTo::SectionTruss::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SectionTruss" << std::endl;
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(SectionBase);
    ar& BOOST_SERIALIZATION_NVP(mArea);
    ar& BOOST_SERIALIZATION_NVP(mAreaParameters);
#ifdef DEBUG_SERIALIZATION
    std::cout << "end serialize SectionTruss" << std::endl;
#endif
}


BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::SectionTruss)
#endif // ENABLE_SERIALIZATION
