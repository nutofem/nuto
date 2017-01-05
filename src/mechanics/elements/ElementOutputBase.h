// $Id $
#pragma once

#include <vector>



namespace NuTo
{
class ElementOutputIpData;


template <class T, int rows, int cols> class FullMatrix;
template <class T, int rows> class FullVector;
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

    virtual FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixDouble();

    virtual BlockFullMatrix<double>& GetBlockFullMatrixDouble();

    virtual FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& GetFullMatrixInt();

    virtual Eigen::VectorXd& GetFullVectorDouble();

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

