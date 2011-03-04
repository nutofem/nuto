// $Id$
#ifndef ELEMENTDATACRACKBASE_H_
#define ELEMENTDATACRACKBASE_H_

#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/cracks/CrackBase.h"

#include <vector>

namespace NuTo
{
//! @author Daniel Arnold
//! @date October 2010
//! @brief class to store the crack informations for the elements
class ElementDataCrackBase : public virtual ElementDataBase
{

public:
    //! @brief constructor
	ElementDataCrackBase();

	virtual ~ElementDataCrackBase();

    //! @brief gets the cracks of an element
    //! @return vector to cracks
	std::vector< CrackBase*>& GetCracks();

    //! @brief gets the number of cracks for an element
    //! @return number of cracks
	int GetNumCracks()const;

	//! @brief adds a crack to the element
	//! @param rCrack  crack
	//! @return the local crack number, the crack is either append to the list, or the existing local number is returned
    unsigned int AddCrack(CrackBase* rCrack);

protected:
    std::vector<CrackBase*> mCracks;
};
}
#endif // ELEMENTDATACRACKBASE_H_
