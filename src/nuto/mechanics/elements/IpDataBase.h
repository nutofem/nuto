// $Id$ 
#ifndef IPDATABASE_H_
#define IPDATABASE_H_

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeBase.h"
#include "nuto/visualize/VisualizeUnstructuredGrid.h"
#include <boost/ptr_container/ptr_list.hpp>
#endif // ENABLE_VISUALIZE

#include <string>
#include <vector>

#include "nuto/math/FullMatrix_Def.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

namespace NuTo
{
class ConstitutiveBase;
class ConstitutiveStaticDataBase;
class ElementBase;
class VisualizeComponentBase;
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

	virtual ConstitutiveStaticDataBase* GetStaticData();

    virtual const ConstitutiveStaticDataBase* GetStaticData()const;

    virtual void SetStaticData(ConstitutiveStaticDataBase* rStaticData);

    virtual void GetLocalIntegrationPointCoordinates2D(boost::array<double,2 >& rLocalCoordinatesFacet)const;

	virtual void SetLocalIntegrationPointCoordinates2D(const boost::array<double,2 >& rLocalCoordinatesFacet);

	virtual void GetLocalIntegrationPointCoordinates3D(boost::array<double,3 >& rLocalCoordinatesFacet)const;

	virtual void SetLocalIntegrationPointCoordinates3D(const boost::array<double,3 >& rLocalCoordinatesFacet);

	virtual double GetIntegrationPointWeight()const;

	virtual void SetIntegrationPointWeight(double);

    //! @brief returns the enum of IP data type
    //! @return enum of IPDataType
    virtual NuTo::IpData::eIpDataType GetIpDataType()const=0;

    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(std::string rFileName, double rMacroLength, double rCoordinates[2], std::string rIPName);

    //! @brief sets the fine scale parameter for all ips
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(const std::string& rName, double rParameter);

    //! @brief sets the fine scale parameter for all ips
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(const std::string& rName, std::string rParameter);

#ifdef ENABLE_VISUALIZE
	//Visualize for all integration points the fine scale structure
	virtual void VisualizeIpMultiscale(VisualizeUnstructuredGrid& rVisualize,
			const boost::ptr_list<NuTo::VisualizeComponentBase>& rWhat, bool rVisualizeDamage)const;
#endif

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

};
}
#endif /* IPDATABASE_H_ */
