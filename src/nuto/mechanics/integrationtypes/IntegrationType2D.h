// $Id$
#ifndef IntegrationType2D_H
#define IntegrationType2D_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
//! @author Jörg F. Unger, ISM
//! @date March 2009
//! @brief ... abstract class for all integration types in 2D
class IntegrationType2D : public IntegrationTypeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType2D() {};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

    //! @brief returns the dimension of the integration type
    //! @return dimension = 1, 2 or 3
    int GetCoordinateDimension()const
    {
        return 2;
    }

    //! @brief ... check compatibility between element type and integration type
    //! @param rElementType ... element type
    //! @return ... <B>true</B> if the element is compatible with the constitutive relationship, <B>false</B> otherwise.
    virtual bool CheckElementCompatibility(NuTo::Element::eElementType rElementType) const;

protected:


};
}

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType2D)
#endif // ENABLE_SERIALIZATION

#endif //IntegrationType2D_H
