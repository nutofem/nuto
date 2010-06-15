// $ld: $ 
#ifndef IPDATAEMPTY_H_
#define IPDATAEMPTY_H_

#include "nuto/mechanics/elements/IpDataBase.h"

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class IpDataEmpty : public IpDataBase
{
public:
	IpDataEmpty() : IpDataBase()
	{}

	void Initialize(const ElementBase* rElement, const ConstitutiveBase* rConstitutive)
	{}

};
}
#endif /* IPDATAEMPTY_H_ */
