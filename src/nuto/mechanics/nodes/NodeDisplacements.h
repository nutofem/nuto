/*
 * NodeDisplacements.h
 *
 *  Created on: 16 Apr 2015
 *      Author: ttitsche
 */

#ifndef NODEDISPLACEMENTS_H_
#define NODEDISPLACEMENTS_H_

#include "nuto/mechanics/nodes/NodeBase.h"



namespace NuTo
{
template<int TNumDisplacements, int TNumTimeDerivatives>
class NodeDisplacements: public virtual NodeBase
{
public:
    NodeDisplacements()
    {
        static_assert(TNumDisplacements >= 0 and TNumDisplacements <= 3, "The node must have 0,1,2 or 3 displacements.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2, "The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mDisplacements[iTimeDerivative] = Eigen::Matrix<double, TNumDisplacements, 1>::Zero();


        mDofDisplacements.fill(-1);
    }

    NodeDisplacements(const NodeDisplacements&) = default;
    NodeDisplacements& operator=(const NodeDisplacements&) = default;

    virtual ~NodeDisplacements()
    {
    }

    int GetNumDisplacements() const override
    {
        return TNumDisplacements;
    }
    int GetDofDisplacement(int rComponent) const override
    {
        assert(rComponent >= 0);
        assert(rComponent < TNumDisplacements);
        return mDofDisplacements[rComponent];
    }
    double GetDisplacement(short rIndex) const override
    {
        assert(rIndex >= 0);
        assert(rIndex < TNumDisplacements);
        return mDisplacements[0][rIndex];
    }

    inline const Eigen::Matrix<double, 1, 1>& GetDisplacements1D(int rTimeDerivative) const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetDisplacements1D] Not implemented for " + std::to_string(TNumDisplacements) + " displacements.");
    }
    inline const Eigen::Matrix<double, 2, 1>& GetDisplacements2D(int rTimeDerivative) const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetDisplacements2D] Not implemented for " + std::to_string(TNumDisplacements) + " displacements.");
    }
    inline const Eigen::Matrix<double, 3, 1>& GetDisplacements3D(int rTimeDerivative) const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetDisplacements3D] Not implemented for " + std::to_string(TNumDisplacements) + " displacements.");
    }
    const Eigen::Matrix<double, Eigen::Dynamic, 1> GetDisplacements(int rTimeDerivative) const override
    {
        assert(rTimeDerivative >= 0);
        assert(rTimeDerivative <= TNumTimeDerivatives);
        return mDisplacements[rTimeDerivative];
    }

    inline void SetDisplacements1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rDisplacements) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetDisplacements1D] Not implemented for " + std::to_string(TNumDisplacements) + " displacements.");
    }
    inline void SetDisplacements2D(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rDisplacements) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetDisplacements2D] Not implemented for " + std::to_string(TNumDisplacements) + " displacements.");
    }
    inline void SetDisplacements3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rDisplacements) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetDisplacements3D] Not implemented for " + std::to_string(TNumDisplacements) + " displacements.");
    }
    void SetDisplacements(int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rDisplacements) override
    {
        assert(rTimeDerivative >= 0);
        assert(rTimeDerivative <= TNumTimeDerivatives);
        mDisplacements[rTimeDerivative] = rDisplacements;
    }

protected:
    std::array<Eigen::Matrix<double, TNumDisplacements, 1>, TNumTimeDerivatives + 1> mDisplacements;
    std::array<int, TNumDisplacements> mDofDisplacements;
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
        ar & boost::serialization::make_array(mDisplacements.data(), mDisplacements.size());
        ar & boost::serialization::make_array(mDofDisplacements.data(), mDofDisplacements.size());
    }
#endif // ENABLE_SERIALIZATION

};

// ****************** GET TNumTimeDerivatives = 0 ************************
template<>
inline const Eigen::Matrix<double, 1, 1>& NodeDisplacements<1,0>::GetDisplacements1D(int rTimeDerivative) const
{
    return mDisplacements[0];
}
template<>
inline const Eigen::Matrix<double, 2, 1>& NodeDisplacements<2,0>::GetDisplacements2D(int rTimeDerivative) const
{
    return mDisplacements[0];
}
template<>
inline const Eigen::Matrix<double, 3, 1>& NodeDisplacements<3,0>::GetDisplacements3D(int rTimeDerivative) const
{
    return mDisplacements[0];
}

// ****************** GET TNumTimeDerivatives = 1 ************************
template<>
inline const Eigen::Matrix<double, 1, 1>& NodeDisplacements<1,1>::GetDisplacements1D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    return mDisplacements[rTimeDerivative];
}
template<>
inline const Eigen::Matrix<double, 2, 1>& NodeDisplacements<2,1>::GetDisplacements2D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    return mDisplacements[rTimeDerivative];
}
template<>
inline const Eigen::Matrix<double, 3, 1>& NodeDisplacements<3,1>::GetDisplacements3D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    return mDisplacements[rTimeDerivative];
}

// ****************** GET TNumTimeDerivatives = 1 ************************
template<>
inline const Eigen::Matrix<double, 1, 1>& NodeDisplacements<1,2>::GetDisplacements1D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    return mDisplacements[rTimeDerivative];
}
template<>
inline const Eigen::Matrix<double, 2, 1>& NodeDisplacements<2,2>::GetDisplacements2D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    return mDisplacements[rTimeDerivative];
}
template<>
inline const Eigen::Matrix<double, 3, 1>& NodeDisplacements<3,2>::GetDisplacements3D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    return mDisplacements[rTimeDerivative];
}




// ****************** SET TNumTimeDerivatives = 0 ************************
template<>
inline void NodeDisplacements<1,0>::SetDisplacements1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rDisplacements)
{
    mDisplacements[0] = rDisplacements;
}template<>
inline void NodeDisplacements<2,0>::SetDisplacements2D(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rDisplacements)
{
    mDisplacements[0] = rDisplacements;
}
template<>
inline void NodeDisplacements<3,0>::SetDisplacements3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rDisplacements)
{
    mDisplacements[0] = rDisplacements;
}

// ****************** SET TNumTimeDerivatives = 1 ************************
template<>
inline void NodeDisplacements<1,1>::SetDisplacements1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rDisplacements)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    mDisplacements[rTimeDerivative] = rDisplacements;
}
template<>
inline void NodeDisplacements<2,1>::SetDisplacements2D(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rDisplacements)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    mDisplacements[rTimeDerivative] = rDisplacements;
}
template<>
inline void NodeDisplacements<3,1>::SetDisplacements3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rDisplacements)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    mDisplacements[rTimeDerivative] = rDisplacements;
}

// ****************** SET TNumTimeDerivatives = 2 ************************
template<>
inline void NodeDisplacements<1,2>::SetDisplacements1D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rDisplacements)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    mDisplacements[rTimeDerivative] = rDisplacements;
}
template<>
inline void NodeDisplacements<2,2>::SetDisplacements2D(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rDisplacements)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    mDisplacements[rTimeDerivative] = rDisplacements;
}
template<>
inline void NodeDisplacements<3,2>::SetDisplacements3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rDisplacements)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    mDisplacements[rTimeDerivative] = rDisplacements;
}


} /* namespace NuTo */


#endif /* NODEDISPLACEMENTS_H_ */
