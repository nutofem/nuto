// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroup.h"


//! @brief constructor
NuTo::ConstraintNodeGroup::ConstraintNodeGroup(const Group<NodeBase>* rGroup) :
        ConstraintBase(), mGroup(rGroup)
{

}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintNodeGroup::GetNumConstraintEquations()const
{
    return mGroup->GetNumMembers();
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintNodeGroup::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroup::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroup::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroup::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroup::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroup::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNodeGroup::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintNodeGroup" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase)
       & BOOST_SERIALIZATION_NVP(mGroup);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintNodeGroup" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintNodeGroup)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintNodeGroup)
#endif // ENABLE_SERIALIZATION
