// $Id $
#ifndef ELEMENT_OUTPUT_FATIGUE_DATA_
#define ELEMENT_OUTPUT_FATIGUE_DATA_

#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/FatigueDataEnum.h"

namespace NuTo
{
//! @author Vitaliy M. Kindrachuk
//! @date Dec 16, 2014
//! @brief ...
class ElementOutputFatigueData : public ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputFatigueData(NuTo::FatigueData::eFatigueDataType rFatigueDataType)
    {
    	mFatigueDataType = rFatigueDataType;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    ElementOutputFatigueData* Clone() const
    {
    	return new ElementOutputFatigueData(*this);
    }

    NuTo::FatigueData::eFatigueDataType GetFatigueDataType()
	{
        return mFatigueDataType;
	}

private:
    NuTo::FatigueData::eFatigueDataType mFatigueDataType;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputFatigueData)
#endif // ENABLE_SERIALIZATION
#endif /* ELEMENT_OUTPUT_FATIGUE_DATA_ */
