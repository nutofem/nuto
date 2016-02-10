/*
 * NodeNonlocalEqPlasticStrain.h
 *
 *  Created on: 27 May 2015
 *      Author: ttitsche
 */

#ifndef NODENONLOCALEQPLASTICSTRAIN_H_
#define NODENONLOCALEQPLASTICSTRAIN_H_

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @brief class for nonlocal equivalent strain dofs
//! @date May 2015
template<int TNumNonlocalEqPlasticStrain, int TNumTimeDerivatives>
class NodeNonlocalEqPlasticStrain: public virtual NodeBase
{
public:
	//! @brief constructor
    NodeNonlocalEqPlasticStrain()
    {
        static_assert(TNumNonlocalEqPlasticStrain == 0 or TNumNonlocalEqPlasticStrain == 2, 	"The node must have 0 or 2 nonlocal equivalent plastic strain.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2, 	"The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mNonlocalEqPlasticStrain[iTimeDerivative] = Eigen::Matrix<double, TNumNonlocalEqPlasticStrain, 1>::Zero();


        mDofNonlocalEqPlasticStrain.fill(-1);
    }

    //! @brief copy constructor
    NodeNonlocalEqPlasticStrain(const NodeNonlocalEqPlasticStrain&) = default;

    //! @brief operator =
    NodeNonlocalEqPlasticStrain& operator=(const NodeNonlocalEqPlasticStrain&) = default;

    //! @brief destructor
    virtual ~NodeNonlocalEqPlasticStrain()
    {
    }

    int GetNumNonlocalEqPlasticStrain() const override
    {
        return TNumNonlocalEqPlasticStrain;
    }

    int GetDofNonlocalEqPlasticStrain(int rComponent) const override
    {
        assert(rComponent >= 0 and rComponent <= 1);
        return mDofNonlocalEqPlasticStrain[rComponent];
    }

    const Eigen::Matrix<double, 2, 1>& GetNonlocalEqPlasticStrain(int rTimeDerivative) const override
    {
        throw MechanicsException("[NuTo::NodeNonlocalEqPlasticStrain::GetNonlocalEqPlasticStrain] Not implemented for " +
                std::to_string(TNumTimeDerivatives) + " time derivatives and " + std::to_string(TNumNonlocalEqPlasticStrain) + " nonlocal eq. plastic strains.");
    }


    void SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain) override
    {
        throw MechanicsException("[NuTo::NodeNonlocalEqPlasticStrain::SetNonlocalEqPlasticStrain] Not implemented for " +
                std::to_string(TNumTimeDerivatives) + " time derivatives and " + std::to_string(TNumNonlocalEqPlasticStrain) + " nonlocal eq. plastic strains.");
    }

protected:
    std::array<Eigen::Matrix<double, TNumNonlocalEqPlasticStrain, 1>, TNumTimeDerivatives + 1> mNonlocalEqPlasticStrain;
    std::array<int, TNumNonlocalEqPlasticStrain> mDofNonlocalEqPlasticStrain;

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
        ar & boost::serialization::make_array(mNonlocalEqPlasticStrain.data(), mNonlocalEqPlasticStrain.size());
        ar & boost::serialization::make_array(mDofNonlocalEqPlasticStrain.data(), mDofNonlocalEqPlasticStrain.size());
    }

#endif // ENABLE_SERIALIZATION

};

template<>
inline const Eigen::Matrix<double, 2, 1>& NodeNonlocalEqPlasticStrain<2,0>::GetNonlocalEqPlasticStrain(int rTimeDerivative) const
{
    return mNonlocalEqPlasticStrain[0];
}

template<>
inline void NodeNonlocalEqPlasticStrain<2,0>::SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain)
{
    mNonlocalEqPlasticStrain[0] = rNonlocalEqPlasticStrain;
}


template<>
inline const Eigen::Matrix<double, 2, 1>& NodeNonlocalEqPlasticStrain<2,1>::GetNonlocalEqPlasticStrain(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    return mNonlocalEqPlasticStrain[rTimeDerivative];
}

template<>
inline void NodeNonlocalEqPlasticStrain<2,1>::SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 1);
    mNonlocalEqPlasticStrain[rTimeDerivative] = rNonlocalEqPlasticStrain;
}

template<>
inline const Eigen::Matrix<double, 2, 1>& NodeNonlocalEqPlasticStrain<2,2>::GetNonlocalEqPlasticStrain(int rTimeDerivative) const
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    return mNonlocalEqPlasticStrain[rTimeDerivative];
}

template<>
inline void NodeNonlocalEqPlasticStrain<2,2>::SetNonlocalEqPlasticStrain(int rTimeDerivative, const Eigen::Matrix<double, 2, 1>& rNonlocalEqPlasticStrain)
{
    assert(rTimeDerivative >= 0);
    assert(rTimeDerivative <= 2);
    mNonlocalEqPlasticStrain[rTimeDerivative] = rNonlocalEqPlasticStrain;
}

} /* namespace NuTo */

#endif /* NODENONLOCALEQPLASTICSTRAIN_H_ */
