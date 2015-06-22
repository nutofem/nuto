// $Id: LoadSurface3D.h 178 2009-12-11 20:53:12Z eckardt4 $
#ifndef LoadSurfaceBase2D_H
#define LoadSurfaceBase2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/loads/LoadBase.h"

namespace NuTo
{
class NodeBase;
class Element2D;
class StructureBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all surface loads in 3D
class LoadSurfaceBase2D : public LoadBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    LoadSurfaceBase2D(int rLoadCase, StructureBase* rStructure, int rElementGroupId, int rNodeGroupId);

    //! @brief adds the load to global sub-vectors
    //! @param rLoadCase number of the current load case
    //! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
    void AddLoadToGlobalSubVectors(int rLoadCase, NuTo::FullVector<double,Eigen::Dynamic>& rActiceDofsLoadVector, NuTo::FullVector<double,Eigen::Dynamic>& rDependentDofsLoadVector)const;

    //! @brief calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates ... global coordinates
    //! @param rNormal ... normal to the surface (pointing outwards)
    //! @param rLoadVector ... load vector
    virtual void CalculateSurfaceLoad(NuTo::FullVector<double,2>& rCoordinates,NuTo::FullVector<double,2>& rNormal,
    		NuTo::FullVector<double,2>& rLoadVector)const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadBase)
           & BOOST_SERIALIZATION_NVP(mElements2D)
           & BOOST_SERIALIZATION_NVP(mIntegrationType2NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType3NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType4NPtr)
           & BOOST_SERIALIZATION_NVP(mIntegrationType5NPtr);
    }
#endif // ENABLE_SERIALIZATION

protected:
    std::vector<std::pair<const Element2D*, int> > mElements2D;
    IntegrationTypeBase* mIntegrationType2NPtr;
    IntegrationTypeBase* mIntegrationType3NPtr;
    IntegrationTypeBase* mIntegrationType4NPtr;
    IntegrationTypeBase* mIntegrationType5NPtr;
};
}//namespace NuTo
#endif //LoadSurfaceBase2D_H

