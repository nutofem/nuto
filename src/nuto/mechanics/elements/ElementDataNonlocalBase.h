// $ld: $ 
#ifndef ELEMENTDATANONLOCALBASE_H_
#define ELEMENTDATANONLOCALBASE_H_

#include "nuto/mechanics/elements/ElementDataBase.h"

#include <vector>

namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 22, 2010
//! @brief class to store the nonlocal weights for regularized constitutive models
class ElementDataNonlocalBase : public virtual ElementDataBase
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
    //! @brief constructor
	ElementDataNonlocalBase();

	virtual ~ElementDataNonlocalBase();

    //! @brief gets the nonlocal elements for a constitutive model
    //! @param rConstitutive constitutive model
    //! @return vector to nonlocal elements
	virtual const std::vector<const ElementBase*>& GetNonlocalElements()const;

    //! @brief gets the number of nonlocal elements for a constitutive model
    //! @param rConstitutive constitutive model
    //! @return number of nonlocal elements
	virtual int GetNumNonlocalElements()const;

    //! @brief adds an element to the nonlocal elements
    //! @param rConstitutive  constitutive model
    //! @return the local element number, the element is either append to the list, or the existing local number is returned
    int AddNonlocalElement(const ElementBase* rElement);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
    std::vector<const ElementBase*> mNonlocalElements;
};
}
#endif /* ELEMENTDATANONLOCALBASE_H_ */
