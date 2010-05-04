// $ld: $ 
#ifndef IPDATABASE_H_
#define IPDATABASE_H_

namespace NuTo
{
class ElementWithDataBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataBase
{
public :

	virtual void Initialize(const ElementWithDataBase* rElement, int rIp)=0;
};
}
#endif /* IPDATABASE_H_ */
