// $Id$
#ifndef LOADNODEGROUPFORCES3D_H
#define LOADNODEGROUPFORCES3D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/loads/LoadNodeGroup.h"

namespace NuTo
{
class NodeBase;
template <class T>
class Group;
template <class T>
class FullMatrix;
;//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... class for all forces applied to a group of nodes in 3D
class LoadNodeGroupForces3D : public LoadNodeGroup
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the force
    //! @param rValue ... value of the force
    LoadNodeGroupForces3D(const Group<NodeBase>* rGroup, const NuTo::FullMatrix<double>& rDirection, double rValue);

    //! @brief adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(NuTo::FullMatrix<double>& rActiceDofsLoadVector, NuTo::FullMatrix<double>& rDependentDofsLoadVector)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadNodeGroup)
        & BOOST_SERIALIZATION_NVP(mValue)
        & BOOST_SERIALIZATION_NVP(mDirection);
    }
#endif // ENABLE_SERIALIZATION

protected:
    double mValue;  //!< prescribed force of the node
    double mDirection[3]; //!< direction of the applied constraint
};
}//namespace NuTo
#endif //LOADNODEGROUPFORCES3D_H

