#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/loads/LoadNode.h"

namespace NuTo
{
class LoadNodeHeatFlux1D : public LoadNode
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief Constructor.
    //! @param rDirection Direction of the heat flux
    //! @param rValue Value of the flux
    LoadNodeHeatFlux1D(int rLoadCase, const NodeBase* rNode, double rDirection, double rValue);

    //! @brief Adds the load to global sub-vectors.
    //! @param rActiceDofsLoadVector Global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector Global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(int rLoadCase,
            NuTo::FullVector<double,Eigen::Dynamic>& rActiceDofsLoadVector,
            NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const;

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
    double mValue;     //!< Prescribed heat flux at the node
    double mDirection; //!< Direction of the flux
};
}
