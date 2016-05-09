// $Id $
#pragma once
#include <boost/assert.hpp>
#include <vector>
#include "nuto/math/FullMatrix_Def.h"
#include "nuto/math/FullVector_Def.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

namespace NuTo
{
class ElementOutputIpData;
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

    virtual BlockFullMatrix<double>& GetBlockFullMatrixDouble();

    virtual FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixInt();

    virtual FullVector<double,Eigen::Dynamic>& GetFullVectorDouble();

    virtual BlockFullVector<double>& GetBlockFullVectorDouble();

    virtual BlockFullVector<int>& GetBlockFullVectorInt();

    virtual FullVector<int,Eigen::Dynamic>& GetFullVectorInt();

    virtual std::vector<int>& GetVectorInt();

    virtual ElementOutputIpData& GetIpData();

    virtual void SetSymmetry(bool rSymmetric);

    virtual bool GetSymmetry()const;

    virtual void SetConstant(bool rConstant);

    virtual bool GetConstant()const;

    virtual ElementOutputBase* Clone() const = 0;
};

}




#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ElementOutputBase)
#endif // ENABLE_SERIALIZATION

