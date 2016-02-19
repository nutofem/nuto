#ifndef NODEWATERVOLUMEFRACTION_H
#define NODEWATERVOLUMEFRACTION_H


#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{

//! @author Volker Hirthammer, BAM
//! @brief class for relative humidity dofs
//! @date May 2015
template<int TNumWaterVolumeFraction, int TNumTimeDerivatives>
class NodeWaterVolumeFraction: public virtual NodeBase
{
public:
    //! @brief constructor
    NodeWaterVolumeFraction()
    {
        static_assert(TNumWaterVolumeFraction >= 0 and TNumWaterVolumeFraction <= 1, 	"The node must have 0 or 1 water volume fraction.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2, 	"The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mWaterVolumeFraction[iTimeDerivative] = Eigen::Matrix<double, TNumWaterVolumeFraction, 1>::Zero();

        mDofWaterVolumeFraction.fill(0);
    }

    //! @brief copy constructor
    NodeWaterVolumeFraction(const NodeWaterVolumeFraction&) = default;

    //! @brief operator =
    NodeWaterVolumeFraction& operator=(const NodeWaterVolumeFraction&) = default;

    //! @brief destructor
    virtual ~NodeWaterVolumeFraction()
    {}

    //! @brief gets the number of water volume dofs
    int GetNumWaterVolumeFraction() const override
    {
        return TNumWaterVolumeFraction;
    }

    //! @brief gets the water volume fraction DOF number
    int GetDofWaterVolumeFraction() const override
    {
        return mDofWaterVolumeFraction[0];
    }

    //! @brief gets the water volume fraction value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    double GetWaterVolumeFraction(int rTimeDerivative) const override
    {
        assert(rTimeDerivative <= 2);
        return mWaterVolumeFraction[rTimeDerivative](0);
    }

    //! @brief sets the water volume fraction value or its time derivatives
    //! @param rTimeDerivative ...time derivative
    //! @param rWaterVolumeFraction ... relative humidity
    void SetWaterVolumeFraction(int rTimeDerivative, double rWaterVolumeFraction) override
    {
        assert(rTimeDerivative <= 2);
        mWaterVolumeFraction[rTimeDerivative](0) = rWaterVolumeFraction;
    }

protected:
    std::array<Eigen::Matrix<double, TNumWaterVolumeFraction, 1>, TNumTimeDerivatives + 1> mWaterVolumeFraction;
    std::array<int, TNumWaterVolumeFraction> mDofWaterVolumeFraction;

private:
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NodeBase);
        ar & boost::serialization::make_array(mWaterVolumeFraction.data(), mWaterVolumeFraction.size());
        ar & boost::serialization::make_array(mDofWaterVolumeFraction.data(), mDofWaterVolumeFraction.size());
    }
#endif // ENABLE_SERIALIZATION

};
} // namespace NuTo


#endif // NODEWATERVOLUMEFRACTION_H
