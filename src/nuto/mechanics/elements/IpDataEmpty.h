// $Id$
#pragma once

#include "nuto/mechanics/elements/IpDataBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class IpDataEmpty : public virtual IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataEmpty();

	void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive);

    //! @brief returns the enum of IP data type
    //! @return enum of IPDataType
    NuTo::IpData::eIpDataType GetIpDataType()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IpDataEmpty)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::IpDataBase, NuTo::IpDataEmpty>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION
