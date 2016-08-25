// $Id $
#pragma once

#include "nuto/math/FullVector_Def.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputFullVectorDouble : public ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputFullVectorDouble()
    {
    };

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    ElementOutputFullVectorDouble* Clone() const
    {
    	return new ElementOutputFullVectorDouble(*this);
    }

    NuTo::FullVector<double,Eigen::Dynamic>& GetFullVectorDouble() override
	{
        return mVector;
	}
private:
    FullVector<double,Eigen::Dynamic> mVector;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputFullVectorDouble)
#endif // ENABLE_SERIALIZATION
