// $Id$
#ifndef ELEMENT_DATA_BASE_H
#define ELEMENT_DATA_BASE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#include <string>
#include <vector>
#include "nuto/mechanics/elements/IpDataEnum.h"

namespace NuTo
{
class ConstitutiveStaticDataBase;
class ConstitutiveBase;
class IntegrationTypeBase;
class ElementBase;

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all possible element data classes
class ElementDataBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ElementDataBase();

    //! @brief constructor
    virtual ~ElementDataBase();

    //! @brief sets the constitutive law for all integration points of the element
    //! @param rConstitutiveLaw constitutive law
    virtual void SetConstitutiveLaw(const ElementBase* rElement, NuTo::ConstitutiveBase* rConstitutiveLaw);

    //! @brief sets the constitutive law for a single integration point of the element
    //! @param rConstitutiveLaw constitutive law
    //! @param rIp integration point
    virtual void SetConstitutiveLaw(const ElementBase* rElement, int rIp, NuTo::ConstitutiveBase* rConstitutiveLaw);

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(int rIp, std::string rFileName);

    //! @brief sets the fine scale parameter for all ips
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(int rIp, const std::string& rName, double rParameter);

    //! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
    //! @param rElement element
    virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement)=0;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    virtual ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    virtual const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    //! @brief returns the constitutive law of an integration point
    //! @param rIp integration point
    //! @return constitutive law
    virtual  ConstitutiveBase* GetConstitutiveLaw(int rIp);

    //! @brief returns the constitutive law of an integration point
    //! @param rIp integration point
    //! @return constitutive law
    virtual  const ConstitutiveBase* GetConstitutiveLaw(int rIp)const;

    //! @brief sets the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @param rElement pointer to element
    //! @param rIntegrationType pointer to integration type
    virtual void SetIntegrationType(const ElementBase* rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

    //! @brief returns a pointer to the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    virtual const IntegrationTypeBase* GetIntegrationType()const;

    //! @brief adds the nonlocal weight to an integration point
    //! @param rLocalIpNumber local Ip
    //! @param rConstitutive constitutive model for which nonlocal data is to be calculated
    //! @param rNonlocalElement element of the nonlocal ip
    //! @param rNonlocalIp local ip number of the nonlocal ip
    //! @param rWeight weight
    virtual void SetNonlocalWeight(int rLocalIpNumber, const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight);

    //! @brief gets the nonlocal elements for a constitutive model
    //! @param rConstitutive constitutive model
    //! @return vector to nonlocal elements
    virtual const std::vector<const NuTo::ElementBase*>& GetNonlocalElements()const;

    //! @brief gets the number of nonlocal elements for a constitutive model
    //! @param rConstitutive constitutive model
    //! @return number of nonlocal elements
    virtual int GetNumNonlocalElements()const;

    //! @brief gets the nonlocal weights
    //! @param rNonlocalElement local element number (should be smaller than GetNonlocalElements().size()
    //! @return vector of weights for all integration points of the nonlocal element
    virtual const std::vector<double>& GetNonlocalWeights(int rIp, int rNonlocalElement)const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

};
}//namespace NuTo

#endif //ELEMENT_DATA_BASE_H

