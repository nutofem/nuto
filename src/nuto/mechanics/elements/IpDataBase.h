#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

#include <string>
#include <vector>

#include <boost/array.hpp>

namespace NuTo
{
class ConstitutiveBase;
class ElementBase;
class VisualizeComponentBase;
class IpDataStaticDataBase;
    namespace Constitutive
    {
        namespace StaticData
        {
            class Component;
        }
    }


namespace IpData
{
    enum class eIpDataType;
}// namespace IpData

#ifdef ENABLE_VISUALIZE
class VisualizeUnstructuredGrid;
#endif // ENABLE_VISUALIZE

//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public :

	virtual ~IpDataBase();

	virtual void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)=0;

	//! @brief adds the weight to an integration point, eventually reallocates the data
	//! @param rNonlocalElement the Element (local number from the nonlocal elements)
	//! @param rNonlocalIp integration point of the nonlocal element
	//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
	//! @param rWeight nonlocal weight
	virtual void SetNonlocalWeight(int rNonlocalElement,int rNonlocalIp, int rNumIps, double rWeight);

	//! @brief delete the nonlocal elements
    //! @param rConstitutive  constitutive model
	virtual void DeleteNonlocalWeights();

	//! @brief return the nonlocal weights
	//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
	//! @return nonlocal weights
	virtual const std::vector<double>& GetNonlocalWeights(int rNonlocalElement)const;

	virtual Constitutive::StaticData::Component* GetConstitutiveStaticData();

    virtual const Constitutive::StaticData::Component* GetConstitutiveStaticData() const;

    virtual void SetConstitutiveStaticData(Constitutive::StaticData::Component* rStaticData);

    virtual IpDataStaticDataBase& GetIpData();

    virtual const IpDataStaticDataBase& GetIpData() const;

    virtual void GetLocalIntegrationPointCoordinates2D(boost::array<double,2 >& rLocalCoordinatesFacet)const;

	virtual void SetLocalIntegrationPointCoordinates2D(const boost::array<double,2 >& rLocalCoordinatesFacet);

	virtual void GetLocalIntegrationPointCoordinates3D(boost::array<double,3 >& rLocalCoordinatesFacet)const;

	virtual void SetLocalIntegrationPointCoordinates3D(const boost::array<double,3 >& rLocalCoordinatesFacet);

    virtual double GetIntegrationPointWeight() const;

    virtual void SetIntegrationPointWeight(double);
    //! @brief returns the enum of IP data type
    //! @return enum of IPDataType
    virtual NuTo::IpData::eIpDataType GetIpDataType()const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

};
}
