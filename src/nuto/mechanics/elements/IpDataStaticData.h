// $ld: $ 
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
public:
	IpDataStaticData();

	virtual ~IpDataStaticData();

	void Initialize(const ElementWithDataBase* rElement, int rIp);

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
#endif /* IPDATASTATICDATA_H_ */
