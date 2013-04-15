// $Id$
#ifndef LOADNODEFORCES3D_H
#define LOADNODEFORCES3D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/loads/LoadNode.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all forces applied to a single node
class LoadNodeForces3D : public LoadNode
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the force
    //! @param rValue ... value of the force
    LoadNodeForces3D(const NodeBase* rNode, const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDirection, double rValue);

    //! @brief adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rActiceDofsLoadVector, NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rDependentDofsLoadVector)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadNode)
        & BOOST_SERIALIZATION_NVP(mValue)
        & BOOST_SERIALIZATION_NVP(mDirection);
    }
#endif // ENABLE_SERIALIZATION

protected:
    double mValue;  //!< prescribed absolute value of the force at the node
    double mDirection[3]; //!< direction of the applied force (normalized)
};
}//namespace NuTo
#endif //LOADNODEFORCES3D_H

