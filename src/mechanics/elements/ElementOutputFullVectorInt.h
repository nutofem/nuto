// $Id $
#pragma once

#include "math/FullVector_Def.h"

#include "mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputFullVectorInt : public ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputFullVectorInt(){};

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    ElementOutputFullVectorInt* Clone() const
    {
    	return new ElementOutputFullVectorInt(*this);
    }

    FullVector<int,Eigen::Dynamic>& GetFullVectorInt() override
	{
        return mVector;
	}

private:
    FullVector<int,Eigen::Dynamic> mVector;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputFullVectorInt)
#endif // ENABLE_SERIALIZATION
