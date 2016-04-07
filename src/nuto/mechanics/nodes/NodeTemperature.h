#pragma once

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{

template<int TNumTemperature, int TNumTimeDerivatives>
class NodeTemperature: public virtual NodeBase
{
public:
    NodeTemperature()
    {
        static_assert(TNumTemperature >= 0 and TNumTemperature <= 1, 	"The node must have 0 or 1 temperatures.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2, 	"The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mTemperature[iTimeDerivative] = Eigen::Matrix<double, TNumTemperature, 1>::Zero();

        mDofTemperature.fill(0);
    }

    //! @brief copy constructor
    NodeTemperature(const NodeTemperature&) = default;

    //! @brief operator =
    NodeTemperature& operator=(const NodeTemperature&) = default;

    //! @brief destructor
    virtual ~NodeTemperature()
    {}

    //! @brief gets the number of temperature dofs
    int GetNumTemperature() const override
    {
        return TNumTemperature;
    }

    //! @brief gets the temperature DOF number
    int GetDofTemperature() const override
    {
        return mDofTemperature[0];
    }

    //! @brief gets the temperature value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    double GetTemperature(int rTimeDerivative) const override
    {
        assert(rTimeDerivative <= 2);
        return mTemperature[rTimeDerivative](0);
    }

    //! @brief sets the temperature value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    //! @param rTemperature ... relative humidity
    void SetTemperature(int rTimeDerivative, double rTemperature) override
    {
        assert(rTimeDerivative <= 2);
        mTemperature[rTimeDerivative](0) = rTemperature;
    }

protected:
    std::array<Eigen::Matrix<double, TNumTemperature, 1>, TNumTimeDerivatives + 1> mTemperature;
    std::array<int, TNumTemperature> mDofTemperature;

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
        ar & boost::serialization::make_array(mTemperature.data(), mTemperature.size());
        ar & boost::serialization::make_array(mDofTemperature.data(), mDofTemperature.size());
    }
#endif // ENABLE_SERIALIZATION

};
}
