// $Id$ 
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
	virtual const std::vector<const ElementBase*>& GetNonlocalElements()const override;

    //! @brief gets the number of nonlocal elements for a constitutive model
    //! @param rConstitutive constitutive model
    //! @return number of nonlocal elements
	virtual int GetNumNonlocalElements()const override;

    //! @brief adds an element to the nonlocal elements
    //! @param rConstitutive  constitutive model
    //! @return the local element number, the element is either append to the list, or the existing local number is returned
    int AddNonlocalElement(const ElementBase* rElement);

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes (saves) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start save ElementDataNonlocalBase" << std::endl;
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase);

        const std::uintptr_t* mNonlocalElementsAdress = reinterpret_cast<const std::uintptr_t*>(mNonlocalElements.data());
        int size = mNonlocalElements.size();
        ar & boost::serialization::make_nvp("mNonlocalElements_size", size);
        ar & boost::serialization::make_nvp("mNonlocalElements", boost::serialization::make_array(mNonlocalElementsAdress, size));
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish save ElementDataNonlocalBase" << std::endl;
    #endif
    }

    //! @brief deserializes(loads) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start load ElementDataNonlocalBase" << std::endl;
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ElementDataBase);

        int size = 0;
        ar & boost::serialization::make_nvp("mNonlocalElements_size", size);
        std::uintptr_t* mNonlocalElementsAdress = new std::uintptr_t[size];
        ar & boost::serialization::make_nvp("mNonlocalElements", boost::serialization::make_array(mNonlocalElementsAdress, size));
        mNonlocalElements.assign(reinterpret_cast<ElementBase**>(&mNonlocalElementsAdress[0]), reinterpret_cast<ElementBase**>(&mNonlocalElementsAdress[size]));
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish load ElementDataNonlocalBase" << std::endl;
    #endif
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    virtual void SetElementPtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mElementMapCast) override;
#endif  // ENABLE_SERIALIZATION

protected:
    std::vector<const ElementBase*> mNonlocalElements;
};
} // namespace NuTo

BOOST_CLASS_EXPORT_KEY(NuTo::ElementDataNonlocalBase)
#endif /* ELEMENTDATANONLOCALBASE_H_ */
