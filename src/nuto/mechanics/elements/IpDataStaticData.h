// $Id$ 
#ifndef IPDATASTATICDATA_H_
#define IPDATASTATICDATA_H_

#include "nuto/mechanics/elements/IpDataStaticDataBase.h"

namespace NuTo
{
class ConstitutiveStaticData;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataStaticData : public IpDataStaticDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataStaticData();

	virtual ~IpDataStaticData();

	void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive);

    //! @brief returns the enum of IP data type
    //! @return enum of IPDataType
    const NuTo::IpData::eIpDataType GetIpDataType()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::IpDataStaticData)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::IpDataStaticDataBase, NuTo::IpDataStaticData>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION
#endif /* IPDATASTATICDATA_H_ */
