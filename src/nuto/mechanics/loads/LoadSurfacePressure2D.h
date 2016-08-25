// $Id: LoadSurface3D.h 178 2009-12-11 20:53:12Z eckardt4 $
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/loads/LoadSurfaceBase2D.h"

namespace NuTo
{
class NodeBase;
class StructureBase;
//! @author Jörg F. Unger, ISM
//! @date August 2013
//! @brief ... class for surface loads in 3D with a const direction and amplitude of the load
class LoadSurfacePressure2D : public LoadSurfaceBase2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    LoadSurfacePressure2D(int rLoadCase, StructureBase* rStructure, int rElementGroupId, int rNodeGroupId,
            double rPressure);

    //! @brief calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates ... global coordinates
    //! @param rNormal ... normal to the surface (pointing outwards)
    //! @param rLoadVector ... load vector
    void CalculateSurfaceLoad(NuTo::FullVector<double,2>& rCoordinates,NuTo::FullVector<double,2>& rNormal,
    		NuTo::FullVector<double,2>& rLoadVector)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadSurfaceBase2D)
           & BOOST_SERIALIZATION_NVP(mPressure);
    }


    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast   std::map containing the old and new Addresses
    virtual void SetElementPtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mElementMapCast) override
    {
            NuTo::LoadSurfaceBase2D::SetElementPtrAfterSerialization(mElementMapCast);
    }
#endif // ENABLE_SERIALIZATION

protected:
    double mPressure;

private:
    //! @brief just for serialization
    LoadSurfacePressure2D(){ }
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LoadSurfacePressure2D)
#endif

