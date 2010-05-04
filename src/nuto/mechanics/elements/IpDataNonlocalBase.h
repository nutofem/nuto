// $ld: $ 
#ifndef IPDATANONLOCALBASE_H_
#define IPDATANONLOCALBASE_H_

#include <vector>
#include "nuto/mechanics/elements/IpDataBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataNonlocalBase : public virtual IpDataBase
{
public:
	IpDataNonlocalBase() : IpDataBase()
	{}

	void SetNonlocalWeights(int rNonlocalElement, const std::vector<double>& rWeights);

	const std::vector<double>& GetNonlocalWeights(int rNonlocalElement)const;

protected:
    //! @brief for each nonlocal element (local number stored at the element level) the nonlocal weights (for each integration point of the nonlocal element)
	std::vector<std::vector<double> > mNonlocalWeights;
};
}
#endif /* IPDATANONLOCALBASE_H_ */
