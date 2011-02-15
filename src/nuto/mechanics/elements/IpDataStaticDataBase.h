// $Id$ 
#ifndef IPDATASTATICDATABASE_H_
#define IPDATASTATICDATABASE_H_

#include "nuto/mechanics/elements/IpDataBase.h"

namespace NuTo
{
class ConstitutiveStaticDataBase;
//! @author Joerg F. Unger
//! @date Apr 28, 2010
//! @brief ...
class IpDataStaticDataBase : public virtual IpDataBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
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
    //! @brief sets the fine scale model (deserialization from a binary file)
    virtual void SetFineScaleModel(std::string rFileName);

    //! @brief sets the fine scale parameter
    //! @parameter rName name of the parameter, e.g. YoungsModulus
    //! @parameter rParameter value of the parameter
    virtual void SetFineScaleParameter(const std::string& rName, double rParameter);

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
