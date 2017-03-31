#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/loads/LoadNodeGroup.h"

namespace NuTo
{
template <class T>
class Group;

//! @brief Class for all forces applied to a group of nodes in 1D
class LoadNodeGroupForces1D : public LoadNodeGroup
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief Constructor
    //! @param direction Direction of the force
    //! @param value Value of the force
    LoadNodeGroupForces1D(const Group<NodeBase>* rGroup, double direction, double value);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(Eigen::VectorXd& rActiceDofsLoadVector,
                                   Eigen::VectorXd& rDependentDofsLoadVector) const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadNodeGroup) & BOOST_SERIALIZATION_NVP(mValue);
    }
#endif // ENABLE_SERIALIZATION

protected:
    double mValue; //!< prescribed load of the node
    double mDirection; //!< direction of the force
};
} // namespace NuTo
