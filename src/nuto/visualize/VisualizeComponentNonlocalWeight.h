// $ld: $ 
#ifndef VISUALIZECOMPONENTNONLOCALWEIGHT_H_
#define VISUALIZECOMPONENTNONLOCALWEIGHT_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/visualize/VisualizeComponentBase.h"

namespace NuTo
{
class ElementBase;

//! @author Joerg F. Unger
//! @date Apr 27, 2010
//! @brief visualize the nonlocal weights for integration point mIp in element mElement
class VisualizeComponentNonlocalWeight : public VisualizeComponentBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
	VisualizeComponentNonlocalWeight(const ElementBase* rElement, int rElementId, int rIp);

    std::string GetComponentName()const;

    int GetElementId()const;

    const NuTo::ElementBase* GetElement()const;

    int GetIp()const;

    inline NuTo::VisualizeBase::eVisualizeWhat GetComponentEnum()const
    {
    	return NuTo::VisualizeBase::NONLOCAL_WEIGHT;
    }

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    VisualizeComponentNonlocalWeight(){};
    int mElementId;
    const ElementBase* mElement;
    int	mIp;
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::VisualizeComponentNonlocalWeight)
#endif // ENABLE_SERIALIZATION

#endif /* VISUALIZECOMPONENTNONLOCALWEIGHT_H_ */
