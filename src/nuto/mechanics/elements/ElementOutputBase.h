// $Id $
#ifndef ELEMENT_OUTPUTBASE_H_
#define ELEMENT_OUTPUTBASE_H_
#include <boost/assert.hpp>
#include <vector>
#include "nuto/math/FullMatrix_Def.h"
#include "nuto/math/FullVector_Def.h"
#include "nuto/mechanics/elements/IpDataEnum.h"


namespace NuTo
{
//! @author Joerg F. Unger
//! @date Apr 29, 2010
//! @brief ...
class ElementOutputBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION
public:
    ElementOutputBase();

    virtual ~ElementOutputBase();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    virtual FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixDouble();

    virtual FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixInt();

    virtual FullVector<double,Eigen::Dynamic>& GetFullVectorDouble();

    virtual FullVector<int,Eigen::Dynamic>& GetFullVectorInt();

    virtual std::vector<int>& GetVectorInt();

    virtual NuTo::IpData::eIpStaticDataType GetIpDataType();

    virtual void SetSymmetry(bool rSymmetric);

    virtual bool GetSymmetry()const;

    virtual void SetConstant(bool rConstant);

    virtual bool GetConstant()const;

    virtual ElementOutputBase* Clone() const = 0;
};

NuTo::ElementOutputBase* new_clone( const ElementOutputBase& o);
}




#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputBase)
#endif // ENABLE_SERIALIZATION
#endif /* ELEMENT_OUTPUTBASE_H_ */
