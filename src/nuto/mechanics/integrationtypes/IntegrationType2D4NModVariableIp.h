// $Id$
#ifndef IntegrationType2D4NModVariableIp_H
#define IntegrationType2D4NModVariableIp_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif //ENABLE_SERIALIZATION
#include <boost/ptr_container/ptr_map.hpp>

#include "nuto/mechanics/integrationtypes/IntegrationPointBase.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2DMod.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NConstVariableIp.h"

namespace NuTo
{
//! @author Daniel Arnold, ISM
//! @date February 2011
//! @brief ... integration types in 2D with four nodes and variable number of integration points
class IntegrationType2D4NModVariableIp : public IntegrationType2DMod
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    IntegrationType2D4NModVariableIp(const std::string rName,const int rNumIp);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:

};
}

#endif //IntegrationType2D4NModVariableIp_H
