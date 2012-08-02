// $Id $
#ifndef ELEMENT_OUTPUT_FULLMATRIX_INT_H_
#define ELEMENT_OUTPUT_FULLMATRIX_INT_H_

#include "nuto/math/FullMatrix.h"

#include "nuto/mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputIpData : public ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputIpData(NuTo::IpData::eIpStaticDataType rIpDataType)
    {
    	mIpDataType = rIpDataType;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    FullMatrix<double>& GetFullMatrixDouble()
	{
        return mMatrix;
	}

    NuTo::IpData::eIpStaticDataType GetIpDataType()
	{
        return mIpDataType;
	}

private:
    NuTo::IpData::eIpStaticDataType mIpDataType;
    FullMatrix<double> mMatrix;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputIpData)
#endif // ENABLE_SERIALIZATION
#endif /* ELEMENT_OUTPUT_FULLMATRIX_INT_H_ */
