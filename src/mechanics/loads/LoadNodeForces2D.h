#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "mechanics/loads/LoadNode.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

namespace NuTo
{

class LoadNodeForces2D : public LoadNode
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief Constructor
    //! @param Irection Direction of the force
    //! @param value Value of the force
    LoadNodeForces2D(const NodeBase* node, const Eigen::MatrixXd& direction, double value);

    //! @brief Adds the load to global sub-vectors
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(StructureOutputBlockVector& externalLoad) const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadNode);
        ar & BOOST_SERIALIZATION_NVP(mValue);
        ar & BOOST_SERIALIZATION_NVP(mDirection);
    }
#endif // ENABLE_SERIALIZATION

protected:
    double mValue;  //!< prescribed absolute value of the force at the node
    double mDirection[2]; //!< direction of the applied force (normalized)

private:
    LoadNodeForces2D() = default;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LoadNodeForces2D)
#endif
