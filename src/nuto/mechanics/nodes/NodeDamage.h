#pragma once

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{

template<int TNumDamage, int TNumTimeDerivatives>
class NodeDamage: public virtual NodeBase
{
public:
    NodeDamage()
    {
        static_assert(TNumDamage >= 0 and TNumDamage <= 1,    "The node must have 0 or 1 damage dofs.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2,    "The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mDamage[iTimeDerivative] = Eigen::Matrix<double, TNumDamage, 1>::Zero();

        mDofDamage.fill(-999);
    }

    //! @brief copy constructor
    NodeDamage(const NodeDamage&) = default;

    //! @brief operator =
    NodeDamage& operator=(const NodeDamage&) = default;

    //! @brief destructor
    virtual ~NodeDamage()
    {}

    //! @brief gets the number of Damage dofs
    int GetNumDamage() const override
    {
        return TNumDamage;
    }

    //! @brief gets the Damage DOF number
    int GetDofDamage() const override
    {
        return mDofDamage[0];
    }

    //! @brief gets the Damage value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    double GetDamage(int rTimeDerivative) const override
    {
        assert(rTimeDerivative <= 2);
        return mDamage[rTimeDerivative](0);
    }

    //! @brief sets the Damage value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    //! @param rDamage ... relative humidity
    void SetDamage(int rTimeDerivative, double rDamage) override
    {
        assert(rTimeDerivative <= 2);
        mDamage[rTimeDerivative](0) = rDamage;
    }

protected:
    std::array<Eigen::Matrix<double, TNumDamage, 1>, TNumTimeDerivatives + 1> mDamage;
    std::array<int, TNumDamage> mDofDamage;

private:
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuTo::NodeBase);
        ar & boost::serialization::make_array(mDamage.data(), mDamage.size());
        ar & boost::serialization::make_array(mDofDamage.data(), mDofDamage.size());
    }
#endif // ENABLE_SERIALIZATION

};
}
