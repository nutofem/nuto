// $Id $
#pragma once


#include "mechanics/elements/ElementOutputBase.h"

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

    ElementOutputFullMatrixInt* Clone() const override
    {
    	return new ElementOutputFullMatrixInt(*this);
    }

    Eigen::MatrixXi& GetFullMatrixInt()override
	{
        return mMatrix;
	}

private:
    Eigen::MatrixXi mMatrix;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputFullMatrixInt)
#endif // ENABLE_SERIALIZATION
