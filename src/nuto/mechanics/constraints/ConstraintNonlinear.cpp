// $Id: ConstraintNonlinear.cpp 328 2010-10-01 14:39:32Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constraints/ConstraintNonlinear.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstraintNonlinear::ConstraintNonlinear()
{
}

// destructor
NuTo::ConstraintNonlinear::~ConstraintNonlinear()
{
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintNonlinear::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinear::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinear::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinear::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinear::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNonlinear::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNonlinear::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintNonlinear" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintNonlinear" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintNonlinear)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintNonlinear)
#endif // ENABLE_SERIALIZATION
