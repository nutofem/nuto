// $Id $
#pragma once

#include "nuto/math/FullMatrix_Def.h"

#include "nuto/mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputFullMatrixInt : public ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputFullMatrixInt(){};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    ElementOutputFullMatrixInt* Clone() const
    {
    	return new ElementOutputFullMatrixInt(*this);
    }

    FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixInt()override
	{
        return mMatrix;
	}

private:
    FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic> mMatrix;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputFullMatrixInt)
#endif // ENABLE_SERIALIZATION
