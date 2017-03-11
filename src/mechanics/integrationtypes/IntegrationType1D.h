#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION

#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
class IntegrationType1D : public IntegrationTypeBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType1D() {}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize IntegrationType1D" << std::endl;
#endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(IntegrationTypeBase);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize IntegrationType1D" << std::endl;
#endif
    }
#endif // ENABLE_SERIALIZATION

    int GetDimension() const override { return 1; }
};
} // namespace nuto

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IntegrationType1D)
#endif


