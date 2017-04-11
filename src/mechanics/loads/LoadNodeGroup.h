#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/loads/LoadBase.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

namespace NuTo
{
template <class T>
class Group;
class NodeBase;

//! @brief Abstract class for all loads applied to a node group
class LoadNodeGroup : public LoadBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    LoadNodeGroup(const Group<NodeBase>* rGroup);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumConstraintEquations() const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadBase) & BOOST_SERIALIZATION_NVP(mGroup);
    }
#endif // ENABLE_SERIALIZATION


protected:
    const Group<NodeBase>* mGroup;
};
} // namespace NuTo
