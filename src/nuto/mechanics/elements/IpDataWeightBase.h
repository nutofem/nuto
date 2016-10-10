// $Id: $
#pragma once

#include "nuto/mechanics/elements/IpDataBase.h"


namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataWeightBase : public virtual IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataWeightBase();

	virtual ~IpDataWeightBase();

	double GetIntegrationPointWeight()const
	{
		return mWeight;
	}

	void SetIntegrationPointWeight(double rWeight)
	{
		mWeight = rWeight ;
	}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
	double mWeight;
};
}
