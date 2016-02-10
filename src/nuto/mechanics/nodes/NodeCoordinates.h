/*
 * NodeCoordinates.h
 *
 *  Created on: 15 Apr 2015
 *      Author: ttitsche
 */

#ifndef NODECOORDINATES_H_
#define NODECOORDINATES_H_

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
template<int TNumCoordinates>
class NodeCoordinates: public virtual NodeBase
{
public:
    NodeCoordinates()
    {
        static_assert(TNumCoordinates >= 0 and TNumCoordinates <= 3, "The node must have 0,1,2 or 3 coordinates.");
        mCoordinates = Eigen::Matrix<double, TNumCoordinates, 1>::Zero();
    }

    NodeCoordinates(const NodeCoordinates&) = default;
    NodeCoordinates& operator=(const NodeCoordinates&) = default;

    virtual ~NodeCoordinates() {}

    int GetNumCoordinates() const override
    {
        return TNumCoordinates;
    }

    double GetCoordinate(short rComponent) const override
    {
        assert(0 <= rComponent and rComponent < TNumCoordinates);
        return mCoordinates[rComponent];
    }
    inline const Eigen::Matrix<double, 1, 1>& GetCoordinates1D() const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetCoordinates1D] Not implemented for " + std::to_string(TNumCoordinates) + " coordinates.");
    }
    inline const Eigen::Matrix<double, 2, 1>& GetCoordinates2D() const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetCoordinates2D] Not implemented for " + std::to_string(TNumCoordinates) + " coordinates.");
    }
    inline const Eigen::Matrix<double, 3, 1>& GetCoordinates3D() const override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::GetCoordinates3D] Not implemented for " + std::to_string(TNumCoordinates) + " coordinates.");
    }
    const Eigen::Matrix<double, Eigen::Dynamic, 1> GetCoordinates() const override
    {
        return mCoordinates;
    }

    inline void SetCoordinates1D(const Eigen::Matrix<double, 1, 1>& rCoordinates) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetCoordinates1D] Not implemented for " + std::to_string(TNumCoordinates) + " coordinates.");
    }
    inline void SetCoordinates2D(const Eigen::Matrix<double, 2, 1>& rCoordinates) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetCoordinates2D] Not implemented for " + std::to_string(TNumCoordinates) + " coordinates.");
    }
    inline void SetCoordinates3D(const Eigen::Matrix<double, 3, 1>& rCoordinates) override
    {
        throw MechanicsException("[NuTo::NodeCoordinates::SetCoordinates3D] Not implemented for " + std::to_string(TNumCoordinates) + " coordinates.");
    }
    void SetCoordinates  (const Eigen::Matrix<double, Eigen::Dynamic, 1>& rCoordinates) override
    {
        assert(TNumCoordinates == rCoordinates.rows());
        mCoordinates = rCoordinates;
    }

protected:
    Eigen::Matrix<double, TNumCoordinates, 1> mCoordinates;

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
        ar & BOOST_SERIALIZATION_NVP(mCoordinates);
    }

#endif // ENABLE_SERIALIZATION

};

//*************************************************
//************  COORDINATES   GET   ***************
//*************************************************

template<>
inline const Eigen::Matrix<double, 1, 1>& NodeCoordinates<1>::GetCoordinates1D() const
{
    return mCoordinates;
}
template<>
inline const Eigen::Matrix<double, 2, 1>& NodeCoordinates<2>::GetCoordinates2D() const
{
    return mCoordinates;
}
template<>
inline const Eigen::Matrix<double, 3, 1>& NodeCoordinates<3>::GetCoordinates3D() const
{
    return mCoordinates;
}

//*************************************************
//************  COORDINATES   SET   ***************
//*************************************************

template<>
inline void NodeCoordinates<1>::SetCoordinates1D(const Eigen::Matrix<double, 1, 1>& rCoordinates)
{
    mCoordinates = rCoordinates;
}
template<>
inline void NodeCoordinates<2>::SetCoordinates2D(const Eigen::Matrix<double, 2, 1>& rCoordinates)
{
    mCoordinates = rCoordinates;
}
template<>
inline void NodeCoordinates<3>::SetCoordinates3D(const Eigen::Matrix<double, 3, 1>& rCoordinates)
{
    mCoordinates = rCoordinates;
}

} /* namespace NuTo */



#endif /* NODECOORDINATES_H_ */
