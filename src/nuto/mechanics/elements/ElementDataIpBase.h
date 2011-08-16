// $Id$ 
#ifndef ELEMENTDATAIPBASE_H_
#define ELEMENTDATAIPBASE_H_

#include <boost/ptr_container/ptr_vector.hpp>

#include "nuto/mechanics/elements/IpDataBase.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

namespace NuTo
{
class IntegrationTypeBase;
class IpDataBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class ElementDataIpBase : public virtual ElementDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    //! @brief constructor
    //! @param rElement			... element for the IP Data
    //! @param rIntegrationType	... integration type
    //! @param rIpDataType		... the IP Data
	ElementDataIpBase(const ElementBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

	virtual ~ElementDataIpBase();

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(int rIp, std::string rFileName, double rLengthCoarseScale, double rCoordinates[2], std::string rIPName);

    //! @brief sets the fine scale parameter for all ips
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(int rIp, const std::string& rName, double rParameter);

    //! @brief sets the fine scale parameter for all ips
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(int rIp, const std::string& rName, std::string rParameter);

#ifdef ENABLE_VISUALIZE
    //Visualize for all integration points the fine scale structure
    void VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
    		const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const;
#endif

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

    //! @brief returns ip data type of the element
    //! implemented with an exception for all element data, reimplementation required for those element data
    //! which actually need an integration type
    //! @return enum to ip data
    NuTo::IpData::eIpDataType GetIpDataType(int  rIp)const;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    //! @brief sets the static data for an integration point of an element
    //! @param rIp integration point
    //! @param rStaticData static data
    virtual void SetStaticData(int rIp, ConstitutiveStaticDataBase* rStaticData);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    //! @brief ...just for serialization
    ElementDataIpBase(){}
    const IntegrationTypeBase *mIntegrationType;
    boost::ptr_vector<IpDataBase> mIpData;
};
}
#endif /* ELEMENTDATAIPBASE_H_ */
