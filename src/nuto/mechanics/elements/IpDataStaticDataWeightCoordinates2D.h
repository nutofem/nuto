// $Id: IpDataStaticDataNonlocal.h 569 2011-08-16 21:13:55Z unger3 $ 
#ifndef IpDataStaticDataWeightCoordinates2D_H_
#define IpDataStaticDataWeightCoordinates2D_H_

#include "nuto/mechanics/elements/IpDataStaticDataBase.h"
#include "nuto/mechanics/elements/IpDataCoordinates2DBase.h"
#include "nuto/mechanics/elements/IpDataWeightBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataStaticDataWeightCoordinates2D : public IpDataStaticDataBase ,public IpDataCoordinates2DBase, public IpDataWeightBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    IpDataStaticDataWeightCoordinates2D();

	virtual ~IpDataStaticDataWeightCoordinates2D();

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
BOOST_CLASS_EXPORT_KEY(NuTo::IpDataStaticDataWeightCoordinates2D)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::IpDataStaticDataBase, NuTo::IpDataStaticDataWeightCoordinates2D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::IpDataCoordinates2DBase, NuTo::IpDataStaticDataWeightCoordinates2D>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::IpDataWeightBase, NuTo::IpDataStaticDataWeightCoordinates2D>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION
#endif /* IpDataStaticDataWeightCoordinates2D_H_ */
