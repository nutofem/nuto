#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/loads/LoadSurfaceBase3D.h"

namespace NuTo
{
class NodeBase;
class StructureBase;

//! @brief Class for surface loads in 3D with a const direction and amplitude of the load
class LoadSurfaceConstDirection3D : public LoadSurfaceBase3D
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief Constructor
    LoadSurfaceConstDirection3D(StructureBase* rStructure, int rElementGroupId, int rNodeGroupId,
                                const Eigen::VectorXd& rLoadVector);

    //! @brief Calculates the surface load as a function of the coordinates and the normal (for pressure)
    //! @param rCoordinates Global coordinates
    //! @param rNormal Normal to the surface (pointing outwards)
    //! @param rLoadVector Load vector
    void CalculateSurfaceLoad(Eigen::Vector3d& rCoordinates, Eigen::Vector3d& rNormal,
                              Eigen::Vector3d& rLoadVector) const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template <class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(LoadSurfaceBase3D) & BOOST_SERIALIZATION_NVP(mLoadVector);
    }
#endif // ENABLE_SERIALIZATION

protected:
    Eigen::Vector3d mLoadVector;
};
} // namespace NuTo
