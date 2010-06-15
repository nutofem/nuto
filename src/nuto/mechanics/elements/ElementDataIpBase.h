// $ld: $ 
#ifndef ELEMENTDATAIPBASE_H_
#define ELEMENTDATAIPBASE_H_

#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#include <boost/ptr_container/serialize_ptr_vector.hpp>

namespace NuTo
{
class IntegrationTypeBase;
class IpDataBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class ElementDataIpBase : public virtual ElementDataBase
{
public:
	ElementDataIpBase(const ElementBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

	virtual ~ElementDataIpBase();

	//! @brief sets the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @param rElement pointer to element
    //! @param rIntegrationType pointer to integration type
    void SetIntegrationType(const ElementBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

    //! @brief returns a pointer to the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    const IntegrationTypeBase* GetIntegrationType()const;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

protected:
    const IntegrationTypeBase *mIntegrationType;
    boost::ptr_vector<IpDataBase> mIpData;
};
}
#endif /* ELEMENTDATAIPBASE_H_ */
