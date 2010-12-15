// $Id: ConstraintLinear.cpp 328 2010-10-01 14:39:32Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constraints/ConstraintLinear.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstraintLinear::ConstraintLinear()
{
}

// destructor
NuTo::ConstraintLinear::~ConstraintLinear()
{
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLinear::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinear::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinear::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinear::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinear::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinear::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinear::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinear" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinear" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinear)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintLinear)
#endif // ENABLE_SERIALIZATION
