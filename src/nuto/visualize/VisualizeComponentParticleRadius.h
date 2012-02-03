// $Id:$
#ifndef VISUALIZECOMPONENTPARTICLERADIUS_H_
#define VISUALIZECOMPONENTPARTICLERADIUS_H_

#include <string>

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"


namespace NuTo
{
//! @author Joerg F. Unger
//! @date Jan, 2012
//! @brief ...
class VisualizeComponentParticleRadius : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    VisualizeComponentParticleRadius();

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::PARTICLE_RADIUS;
    }

    inline std::string GetComponentName()const
    {
    	return std::string("ParticleRadius");
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentParticleRadius)
#endif // ENABLE_SERIALIZATION
#endif /* VISUALIZECOMPONENTPARTICLERADIUS_H_ */
