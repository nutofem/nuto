// $Id$ 
#ifndef IPDATASTATICDATANONLOCAL_H_
#define IPDATASTATICDATANONLOCAL_H_

#include "nuto/mechanics/elements/IpDataStaticDataBase.h"
#include "nuto/mechanics/elements/IpDataNonlocalBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataStaticDataNonlocal : public IpDataStaticDataBase ,public IpDataNonlocalBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataStaticDataNonlocal();

	virtual ~IpDataStaticDataNonlocal();

	virtual void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive);

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
BOOST_CLASS_EXPORT_KEY(NuTo::IpDataStaticDataNonlocal)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::IpDataStaticDataBase, NuTo::IpDataStaticDataNonlocal>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::IpDataNonlocalBase, NuTo::IpDataStaticDataNonlocal>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION
#endif /* IPDATASTATICDATANONLOCAL_H_ */
