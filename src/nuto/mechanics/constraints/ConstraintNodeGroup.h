// $Id$

#ifndef CONSTRAINTNODEGROUP_H
#define CONSTRAINTNODEGROUP_H

#include "nuto/mechanics/constraints/ConstraintEnum.h"

namespace NuTo
{
template <class T>
class Group;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all constraints applied to a node group
class ConstraintNodeGroup
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintNodeGroup(const Group<NodeBase>* rGroup);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION


protected:
    //! @brief just for serialization
    ConstraintNodeGroup(){};
    const Group<NodeBase>* mGroup;
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
//BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintNodeGroup)
#endif // ENABLE_SERIALIZATION
#endif //CONSTRAINTNODEGROUP_H

