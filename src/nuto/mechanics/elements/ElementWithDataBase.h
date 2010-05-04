// $Id: $
#ifndef ElementWithDataBase_H
#define ElementWithDataBase_H

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

namespace NuTo
{
class ConstitutiveBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all finite elements with element data
//! like materials, static data etc.
//! This class has been introduced in order to allow elements without any static data, whose
//! material law is defined e.g. at the structural level
class ElementWithDataBase : public ElementBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:

    //! @brief constructor
	ElementWithDataBase(const StructureBase* rStructure, ElementData::eElementDataType rElementDataType,
			IntegrationType::eIntegrationType rIntegrationType, IpData::eIpDataType rIpDataType);

	//! @brief destructor
	virtual ~ElementWithDataBase();

	//! @brief returns a pointer to the constitutive law for an integration point
    //! @param integration point number (counting from zero)
    //! @return pointer to constitutive law
    const ConstitutiveBase* GetConstitutiveLaw(int rIp)const;

    //! @brief returns a pointer to the constitutive law for an integration point
    //! @param integration point number (counting from zero)
    //! @return pointer to constitutive law
    ConstitutiveBase* GetConstitutiveLaw(int rIp);

    //! @brief sets the constitutive law for an element
    //! @param rConstitutiveLaw Pointer to constitutive law entry
    void SetConstitutiveLaw(ConstitutiveBase* rConstitutiveLaw);

    //! @brief sets the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @param rIntegrationType pointer to integration type
    void SetIntegrationType(const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

    //! @brief returns a pointer to the integration type of an element
    //! implemented with an exception for all elements, reimplementation required for those elements
    //! which actually need an integration type
    //! @return pointer to integration type
    const IntegrationTypeBase* GetIntegrationType()const;

    //! @brief returns the number of integration points
    //! @return number of integration points
    int GetNumIntegrationPoints()const;

    //! @brief returns the weight of an integration point
    //! @param rIpNum integration point
    //! @return weight
    double GetIntegrationPointWeight(int rIpNum)const;

    //! @brief returns the coordinates of an integration point
    //! @param rIpNum integration point
    //! @param rCoordinates coordinates to be returned
    virtual void GetGlobalIntegrationPointCoordinates(int rIpNum, double rCoordinates[3])const=0;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data pointer
    NuTo::ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data pointer
    NuTo::ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief returns the enum of the standard integration type for this element
    virtual NuTo::IntegrationType::eIntegrationType GetStandardIntegrationType()=0;

    //! @brief adds the nonlocal weight to an integration point
    //! @param rLocalIpNumber local Ip
    //! @param rConstitutive constitutive model for which nonlocal data is to be calculated
    //! @param rNonlocalElement element of the nonlocal ip
    //! @param rNonlocalIp local ip number of the nonlocal ip
    //! @param rWeight weight
    void AddNonlocalIp(int rLocalIpNumber, const ConstitutiveBase* rConstitutive,
    		const ElementWithDataBase* rNonlocalElement, int rNonlocalIp, double rWeight);

    //! @brief calculates the volume of an integration point (weight * detJac)
    //! @param rVolume  vector for storage of the ip volumes (area in 2D, length in 1D)
    virtual void GetIntegrationPointVolume(std::vector<double>& rVolume)const=0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementBase)
           & BOOST_SERIALIZATION_NVP(mElementData);
    }
#endif  // ENABLE_SERIALIZATION
#ifdef ENABLE_VISUALIZE
    void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<NuTo::VisualizeComponentBase*>& rWhat) const;
#endif // ENABLE_VISUALIZE

protected:
    ElementDataBase *mElementData;
};

} // namespace NuTo

#endif //ElementWithDataBase_H
