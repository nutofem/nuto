/*
 * NodeNonlocalEqStrain.h
 *
 *  Created on: 25 May 2015
 *      Author: phuschke
 */

#ifndef NODENONLOCALEQSTRAIN_H_
#define NODENONLOCALEQSTRAIN_H_

#include "nuto/mechanics/nodes/NodeBase.h"

namespace NuTo
{
//! @brief class for nonlocal equivalent strain dofs
//! @date May 2015
template<int TNumNonlocalEqStrain, int TNumTimeDerivatives>
class NodeNonlocalEqStrain: public virtual NodeBase
{
public:
	//! @brief constructor
    NodeNonlocalEqStrain()
    {
        static_assert(TNumNonlocalEqStrain >= 0 and TNumNonlocalEqStrain <= 1, 	"The node must have 0 or 1 nonlocal equivalent strain.");
        static_assert(TNumTimeDerivatives >= 0 and TNumTimeDerivatives <= 2, 	"The node must have 0,1 or 2 time derivatives.");
        for (int iTimeDerivative = 0; iTimeDerivative < TNumTimeDerivatives + 1; ++iTimeDerivative)
            mNonlocalEqStrain[iTimeDerivative] = Eigen::Matrix<double, TNumNonlocalEqStrain, 1>::Zero();


        mDofNonlocalEqStrain.fill(-999);
    }

    //! @brief copy constructor
    NodeNonlocalEqStrain(const NodeNonlocalEqStrain&) = default;

    //! @brief operator =
    NodeNonlocalEqStrain& operator=(const NodeNonlocalEqStrain&) = default;

    //! @brief destructor
    virtual ~NodeNonlocalEqStrain()
    {
    }

    int GetNumNonlocalEqStrain() const override
    {
        return TNumNonlocalEqStrain;
    }

    int GetDofNonlocalEqStrain() const override
    {
        return mDofNonlocalEqStrain[0];
    }

    double GetNonlocalEqStrain(int rTimeDerivative) const override
    {
    	assert(rTimeDerivative <= 2);
        return mNonlocalEqStrain[rTimeDerivative](0);
    }


    void SetNonlocalEqStrain(int rTimeDerivative, double rNonlocalEqStrain) override
    {
    	assert(rTimeDerivative <= 2);
        mNonlocalEqStrain[rTimeDerivative](0) = rNonlocalEqStrain;
    }

protected:
    std::array<Eigen::Matrix<double, TNumNonlocalEqStrain, 1>, TNumTimeDerivatives + 1> mNonlocalEqStrain;
    std::array<int, TNumNonlocalEqStrain> mDofNonlocalEqStrain;

private:
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

};


} /* namespace NuTo */
#endif /* NODENONLOCALEQSTRAIN_H_ */
