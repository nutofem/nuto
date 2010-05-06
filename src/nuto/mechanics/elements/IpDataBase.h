// $ld: $ 
#ifndef IPDATABASE_H_
#define IPDATABASE_H_

#include <vector>

namespace NuTo
{
class ConstitutiveBase;
class ElementWithDataBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataBase
{
public :

	virtual ~IpDataBase();

	virtual void Initialize(const ElementWithDataBase* rElement, const ConstitutiveBase* rConstitutive)=0;

	//! @brief adds the weight to an integration point, eventually reallocates the data
	//! @param rNonlocalElement the Element (local number from the nonlocal elements)
	//! @param rNonlocalIp integration point of the nonlocal element
	//! @param rNumIps number of integration points of the nonlocal element (for allocation purpose of not existing)
	//! @param rWeight nonlocal weight
	virtual void SetNonlocalWeight(int rNonlocalElement,int rNonlocalIp, int rNumIps, double rWeight);

	//! @brief return the nonlocal weights
	//! @param rNonlocalElement nonlocal element (between 0 and nonlocal elements.size stored in nonlocal element data)
	//! @return nonlocal weights
	virtual const std::vector<double>& GetNonlocalWeights(int rNonlocalElement)const;
};
}
#endif /* IPDATABASE_H_ */
