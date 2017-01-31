// $Id $
#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>

namespace NuTo
{
class ElementOutputIpData;


template <typename T> class BlockFullMatrix;
template <typename T> class BlockFullVector;

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

    virtual Eigen::MatrixXd& GetFullMatrixDouble();

    virtual BlockFullMatrix<double>& GetBlockFullMatrixDouble();

    virtual Eigen::MatrixXi& GetFullMatrixInt();

    virtual Eigen::VectorXd& GetFullVectorDouble();

    virtual BlockFullVector<double>& GetBlockFullVectorDouble();

    virtual BlockFullVector<int>& GetBlockFullVectorInt();

    virtual Eigen::VectorXi& GetFullVectorInt();

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
