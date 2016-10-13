// $Id$
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#endif  // ENABLE_SERIALIZATION


#include <string>
#include <vector>
#include <boost/array.hpp>



namespace NuTo
{

class ConstitutiveStaticDataBase;
class ConstitutiveBase;
class CrackBase;
class IntegrationTypeBase;
class InterpolationType;
class ElementBase;
class IpDataBase;
class VisualizeComponentBase;
class IpDataStaticDataBase;

namespace ElementData
{
    enum class eElementDataType;
}// namespace ElementData

namespace IpData
{
    enum class eIpDataType;
}// namespace IpData


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

    //! @brief returns true, if the constitutive law has been assigned
    virtual bool HasConstitutiveLawAssigned(int rIp)const;

    //! @brief returns the number of integration points
    //! @return number of integration points
    virtual int GetNumIntegrationPoints()const;

    //! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
    //! @param rElement element
    virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement)=0;

     //! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
     //! @param rElement element
     //! @param rIp Ip
    virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement,int rIp)=0;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    virtual ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    virtual const ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    virtual IpDataStaticDataBase& GetStaticDataBase(int rIp);

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data
    virtual const IpDataStaticDataBase& GetStaticDataBase(int rIp)const;


    //! @brief sets the static data for an integration point of an element
    //! @param rIp integration point
    //! @param rStaticData static data
    virtual void SetStaticData(int rIp, ConstitutiveStaticDataBase* rStaticData);

    //! @brief returns the local coordinate of an integration point
    //! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
    //! there is no standard integration type for the element and the local coordinates are stored directly at the ip
    //! @param rIpNum number of the integration point
    //! @return rLocalCoordinates
    virtual void GetLocalIntegrationPointCoordinates2D(int rIpNum, boost::array<double,2 >& rLocalCoordinates)const;

    //! @brief sets the local coordinate of an integration point
    //! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
    //! there is no standard integration type for the element and the local coordinates are stored directly at the ip
    //! @param rIpNum number of the integration point
    //! @param rLocalCoordinates
    virtual void SetLocalIntegrationPointCoordinates2D(int rIpNum, const boost::array<double,2 >& rLocalCoordinates);

    //! @brief returns the local coordinate of an integration point
    //! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
    //! there is no standard integration type for the element and the local coordinates are stored directly at the ip
    //! @param rIpNum number of the integration point
    //! @return rLocalCoordinates
    virtual void GetLocalIntegrationPointCoordinates3D(int rIpNum, boost::array<double,3 >& rLocalCoordinates)const;

    //! @brief sets the local coordinate of an integration point
    //! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
    //! there is no standard integration type for the element and the local coordinates are stored directly at the ip
    //! @param rIpNum number of the integration point
    //! @return rLocalCoordinates
    virtual void SetLocalIntegrationPointCoordinates3D(int rIpNum, const boost::array<double,3 >& rLocalCoordinates);

    //! @brief gets the weight of an integration point
    //! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
    //! there is no standard integration type for the element and the local coordinates are stored directly at the ip
    //! @param rIpNum number of the integration point
    //! @return weight
    virtual double GetIntegrationPointWeight(int rIpNum)const;

    //! @brief sets the local coordinate of an integration point
    //! usually, it is easier to use the integration type, but for some problems (lattice, XFEM)
    //! there is no standard integration type for the element and the local coordinates are stored directly at the ip
    //! @param rIpNum number of the integration point
    //! @param weight
	virtual void SetIntegrationPointWeight(int rIpNum, double rWeight);

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



    //! @brief returns ip data type of the element
    //! implemented with an exception for all element data, reimplementation required for those element data
    //! which actually need an integration type
    //! @return enum to ip data
    virtual NuTo::IpData::eIpDataType GetIpDataType(int  rIp)const;

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

    //! @brief delete the nonlocal elements
    virtual void DeleteNonlocalElements();

    //! @brief gets the number of nonlocal elements for a constitutive model
    //! @param rConstitutive constitutive model
    //! @return number of nonlocal elements
    virtual int GetNumNonlocalElements()const;

    //! @brief gets the nonlocal weights
    //! @param rNonlocalElement local element number (should be smaller than GetNonlocalElements().size()
    //! @return vector of weights for all integration points of the nonlocal element
    virtual const std::vector<double>& GetNonlocalWeights(int rIp, int rNonlocalElement)const;

    //! @brief returns the enum of element data type
    //! @return enum of ElementDataType
    virtual const NuTo::ElementData::eElementDataType GetElementDataType()const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    //! @brief ElementBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the ElementBase-Pointer is done by searching and casting back the Address in the map
    //! @param mElementMapCast   std::map containing the old and new Addresses
    virtual void SetElementPtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mElementMapCast)
    {
        (void)mElementMapCast;
        /* Do nothing until needed, see e.g. LoadSurfaceBase2D-class*/
    }
#endif  // ENABLE_SERIALIZATION


protected:

};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementDataBase)
#endif // ENABLE_SERIALIZATION

