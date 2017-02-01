#pragma once

#include <functional>

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "mechanics/loads/LoadSurfaceBase2D.h"

namespace NuTo
{
class NodeBase;
class StructureBase;
//! @author Peter Otto, BAM
//! @date August 2016
//! @brief ... class for surface loads in 2D with a function as load
class LoadSurfacePressureFunction2D : public LoadSurfaceBase2D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    LoadSurfacePressureFunction2D(int rLoadCase,
                                  StructureBase* rStructure,
                                  int rElementGroupId,
                                  int rNodeGroupId,
                                  const std::function<Eigen::Vector2d(Eigen::Vector2d)>& rLoadFunction);

    //! @brief calculates the surface load as a function of the coordinates
    //! @param rCoordinates ... global coordinates
    //! @param rNormal ... normal to the surface (pointing outwards)
    //! @param rLoadVector ... load vector
    void CalculateSurfaceLoad(Eigen::Vector2d& rCoordinates,
                              Eigen::Vector2d& rNormal,
                              Eigen::Vector2d& rLoadVector) const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadSurfaceBase2D)
           & BOOST_SERIALIZATION_NVP(mLoadFunction);
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
    std::function<Eigen::Vector2d(Eigen::Vector2d)> mLoadFunction;

private:
    //! @brief just for serialization
    LoadSurfacePressureFunction2D(){ }
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LoadSurfacePressureFunction2D)
#endif
