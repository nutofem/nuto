// $ld: $ 
#ifndef IPDATASTATICDATABASE_H_
#define IPDATASTATICDATABASE_H_

#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/elements/IpDataBase.h"

namespace NuTo
{
class ConstitutiveStaticDataBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataStaticDataBase : public virtual IpDataBase
{
public:
	IpDataStaticDataBase();

	virtual ~IpDataStaticDataBase();

    inline ConstitutiveStaticDataBase* GetStaticData()
	{
        return mStaticData;
	}

    inline const ConstitutiveStaticDataBase* GetStaticData()const
	{
        return mStaticData;
	}

    inline void SetStaticData(ConstitutiveStaticDataBase* rStaticData)
	{
        mStaticData = rStaticData;
	}

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:
	ConstitutiveStaticDataBase* mStaticData;
};
}
#endif /* IPDATASTATICDATABASE_H_ */
