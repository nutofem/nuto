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
	ElementWithDataBase(const StructureBase* rStructure, ElementDataBase::eElementDataType rElementDataType, IntegrationTypeBase::eIntegrationType rIntegrationType);

	//! @brief destructor
	~ElementWithDataBase();

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
    void SetIntegrationType(const NuTo::IntegrationTypeBase* rIntegrationType);

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

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data pointer
    NuTo::ConstitutiveStaticDataBase* GetStaticData(int rIp)const;

    //! @brief returns the static data of an integration point
    //! @param rIp integration point
    //! @return static data pointer
    NuTo::ConstitutiveStaticDataBase* GetStaticData(int rIp);

    //! @brief returns the enum of the standard integration type for this element
    virtual NuTo::IntegrationTypeBase::eIntegrationType GetStandardIntegrationType()=0;

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
    void Visualize(VisualizeUnstructuredGrid& rVisualize, const std::map<std::string,NuTo::VisualizeBase::eVisualizeWhat>& rWhat) const;
#endif // ENABLE_VISUALIZE

protected:
    ElementDataBase *mElementData;
};

} // namespace NuTo

#endif //ElementWithDataBase_H
