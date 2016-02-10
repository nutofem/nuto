/*
 * NodeRotations.h
 *
 *  Created on: 20 Apr 2015
 *      Author: ttitsche
 */

#ifndef NODEROTATIONS_H_
#define NODEROTATIONS_H_

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{

template<int TNumRotations, int TNumTimeDerivatives>
class NodeRotations: public virtual NodeBase
{
public:
    NodeRotations()
    {
        static_assert(TNumRotations == 0 or TNumRotations == 1 or TNumRotations == 3, "The node must have 0,1 or 3 rotations.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2, "The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mRotations[iTimeDerivative] = Eigen::Matrix<double, TNumRotations, 1>::Zero();

        mDofRotations.fill(-1);
    }
    virtual ~NodeRotations() {};

    int GetNumRotations() const override
    {
        return TNumRotations;
    }
    int GetDofRotation(int rComponent) const override
    {
        assert(rComponent >= 0);
        assert(rComponent < TNumRotations);
        return mDofRotations[rComponent];
    }
    double GetRotation(short rIndex) const override
    {
        assert(rIndex >= 0);
        assert(rIndex < TNumRotations);
        return mRotations[0][rIndex];
    }

    inline const Eigen::Matrix<double, 1, 1>& GetRotations2D(int rTimeDerivative) const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetRotations2D] Not implemented for " + std::to_string(TNumRotations) + " rotations.");
    }
    inline const Eigen::Matrix<double, 3, 1>& GetRotations3D(int rTimeDerivative) const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetRotations3D] Not implemented for " + std::to_string(TNumRotations) + " rotations.");
    }
    const Eigen::Matrix<double, Eigen::Dynamic, 1> GetRotations(int rTimeDerivative) const override
    {
        assert(rTimeDerivative >= 0);
        assert(rTimeDerivative <= TNumTimeDerivatives);
        return mRotations[rTimeDerivative];
    }

    inline void SetRotations2D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rRotations) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetRotations2D] Not implemented for " + std::to_string(TNumRotations) + " rotations.");
    }
    inline void SetRotations3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rRotations) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetRotations3D] Not implemented for " + std::to_string(TNumRotations) + " rotations.");
    }
    void SetRotations(int rTimeDerivative, const Eigen::Matrix<double, Eigen::Dynamic, 1>& rRotations) override
    {
        assert(rTimeDerivative >= 0);
        assert(rTimeDerivative <= TNumTimeDerivatives);
        mRotations[rTimeDerivative] = rRotations;
    }

protected:
    std::array<Eigen::Matrix<double, TNumRotations, 1>, TNumTimeDerivatives + 1> mRotations;
    std::array<int, TNumRotations> mDofRotations;

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
        ar & boost::serialization::make_array(mRotations.data(), mRotations.size());
        ar & boost::serialization::make_array(mDofRotations.data(), mDofRotations.size());
    }
#endif // ENABLE_SERIALIZATION
};



// ****************** GET TNumTimeDerivatives = 0 ************************
template<>
inline const Eigen::Matrix<double, 1, 1>& NodeRotations<1,0>::GetRotations2D(int rTimeDerivative) const
{
    return mRotations[0];
}
template<>
inline const Eigen::Matrix<double, 3, 1>& NodeRotations<3,0>::GetRotations3D(int rTimeDerivative) const
{
    return mRotations[0];
}

// ****************** GET TNumTimeDerivatives = 1 ************************
template<>
inline const Eigen::Matrix<double, 1, 1>& NodeRotations<1,1>::GetRotations2D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    return mRotations[rTimeDerivative];
}
template<>
inline const Eigen::Matrix<double, 3, 1>& NodeRotations<3,1>::GetRotations3D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    return mRotations[rTimeDerivative];
}

// ****************** GET TNumTimeDerivatives = 2 ************************
template<>
inline const Eigen::Matrix<double, 1, 1>& NodeRotations<1,2>::GetRotations2D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    return mRotations[rTimeDerivative];
}
template<>
inline const Eigen::Matrix<double, 3, 1>& NodeRotations<3,2>::GetRotations3D(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    return mRotations[rTimeDerivative];
}




// ****************** SET TNumTimeDerivatives = 0 ************************
template<>
inline void NodeRotations<1,0>::SetRotations2D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rRotations)
{
    mRotations[0] = rRotations;
}
template<>
inline void NodeRotations<3,0>::SetRotations3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rRotations)
{
    mRotations[0] = rRotations;
}

// ****************** SET TNumTimeDerivatives = 1 ************************
template<>
inline void NodeRotations<1,1>::SetRotations2D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rRotations)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    mRotations[rTimeDerivative] = rRotations;
}
template<>
inline void NodeRotations<3,1>::SetRotations3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rRotations)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    mRotations[rTimeDerivative] = rRotations;
}

// ****************** SET TNumTimeDerivatives = 2 ************************
template<>
inline void NodeRotations<1,2>::SetRotations2D(int rTimeDerivative, const Eigen::Matrix<double, 1, 1>& rRotations)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    mRotations[rTimeDerivative] = rRotations;
}
template<>
inline void NodeRotations<3,2>::SetRotations3D(int rTimeDerivative, const Eigen::Matrix<double, 3, 1>& rRotations)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    mRotations[rTimeDerivative] = rRotations;
}


} /* namespace NuTo */


#endif /* NODEROTATIONS_H_ */
