// $ld: $ 
#ifndef ELEMENTDATACONSTITUTIVEIPNONLOCAL_H_
#define ELEMENTDATACONSTITUTIVEIPNONLOCAL_H_

#include "nuto/mechanics/elements/ElementDataConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementDataNonlocalBase.h"
#include "nuto/mechanics/elements/ElementDataIpBase.h"

namespace NuTo
{

//! @author Joerg F. Unger
//! @date Apr 23, 2010
//! @brief ...
class ElementDataConstitutiveIpNonlocal : public ElementDataConstitutiveBase, public ElementDataNonlocalBase, public ElementDataIpBase
{

#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
	//! @brief constructor
	ElementDataConstitutiveIpNonlocal(const ElementBase *rElement, const NuTo::IntegrationTypeBase* rIntegrationType, NuTo::IpData::eIpDataType rIpDataType);

	virtual ~ElementDataConstitutiveIpNonlocal();

	//! @brief updates the data related to changes of the constitutive model (e.g. reallocation of static data, nonlocal weights etc.)
    //! @param rElement element
    virtual void InitializeUpdatedConstitutiveLaw(const ElementBase* rElement);

    //! @brief adds the nonlocal weight to an integration point
    //! @param rLocalIpNumber local Ip
    //! @param rConstitutive constitutive model for which nonlocal data is to be calculated
    //! @param rNonlocalElement element of the nonlocal ip
    //! @param rNonlocalIp local ip number of the nonlocal ip
    //! @param rWeight weight
     void SetNonlocalWeight(int rLocalIpNumber,const ElementBase* rNonlocalElement, int rNonlocalIp, double rWeight);

     //! @brief gets the nonlocal weights
     //! @param rNonlocalElement local element number (should be smaller than GetNonlocalElements().size()
     //! @return vector of weights for all integration points of the nonlocal element
     const std::vector<double>& GetNonlocalWeights(int rIp, int rNonlocalElement)const;

#ifdef ENABLE_SERIALIZATION
	//! @brief serializes the class
	//! @param ar         archive
	//! @param version    version
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
	//! @brief ... just for serialization
	ElementDataConstitutiveIpNonlocal(){};
};
}
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementDataConstitutiveIpNonlocal)
namespace boost{
template<>
struct is_virtual_base_of<NuTo::ElementDataConstitutiveBase, NuTo::ElementDataConstitutiveIpNonlocal>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::ElementDataNonlocalBase, NuTo::ElementDataConstitutiveIpNonlocal>: public mpl::true_ {};
template<>
struct is_virtual_base_of<NuTo::ElementDataIpBase, NuTo::ElementDataConstitutiveIpNonlocal>: public mpl::true_ {};
}
#endif // ENABLE_SERIALIZATION

#endif /* ELEMENTDATACONSTITUTIVEIPNONLOCAL_H_ */
