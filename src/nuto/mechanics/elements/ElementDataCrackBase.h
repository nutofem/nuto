// $Id$
#pragma once

#include "nuto/mechanics/elements/ElementDataBase.h"


#include <vector>

namespace NuTo
{
class CrackBase;

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

    //! @brief Set the information that the element is already cracked or not
    //! @param bool (Input) cracked or not
    void IsCracked(const bool rIsCracked);

    //! @brief Give the information if the element is already cracked or not
    //! @return bool cracked or not
    const bool IsCracked() const;

	//! @brief adds a crack to the element
	//! @param rCrack  crack
	//! @return the local crack number, the crack is either append to the list, or the existing local number is returned
    unsigned int AddCrack(CrackBase* rCrack);

protected:
    std::vector<CrackBase*> mCracks;
    bool isCracked;
};
}
