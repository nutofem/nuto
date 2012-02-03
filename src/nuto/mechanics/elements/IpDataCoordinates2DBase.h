// $Id: $
#ifndef IpDataCoordinates2DBase_H_
#define IpDataCoordinates2DBase_H_

#include "nuto/mechanics/elements/IpDataBase.h"
#include "boost/array.hpp"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataCoordinates2DBase : public virtual IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	IpDataCoordinates2DBase();

	virtual ~IpDataCoordinates2DBase();

	void GetLocalIntegrationPointCoordinates2D(boost::array<double,2 >& rLocalCoordinatesFacet)const
	{
		rLocalCoordinatesFacet = mLocalCoordinates;
	}

	void SetLocalIntegrationPointCoordinates2D(const boost::array<double,2 >& rLocalCoordinatesFacet)
	{
		mLocalCoordinates = rLocalCoordinatesFacet ;
	}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
	boost::array<double,2> mLocalCoordinates;
};
}
#endif /* IpDataCoordinates2DBase_H_ */
