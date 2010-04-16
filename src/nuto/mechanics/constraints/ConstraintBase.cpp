// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/constraints/ConstraintBase.h"

//! @brief constructor
NuTo::ConstraintBase::ConstraintBase()
{
}

// destructor
NuTo::ConstraintBase::~ConstraintBase()
{
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintBase::serialize(Archive & ar, const unsigned int version)
{
}
#endif // ENABLE_SERIALIZATION
