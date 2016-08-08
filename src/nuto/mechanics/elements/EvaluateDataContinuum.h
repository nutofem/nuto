#pragma once

#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"

namespace NuTo
{


template<int TDim>
struct EvaluateDataContinuum
{

public:
    // Nodal Values
    // --------------------------------------------------------------------------------------------

    std::map<Node::eDof, Eigen::VectorXd> mNodalValues;
    std::map<Node::eDof, Eigen::VectorXd> mNodalValues_dt1;


    // Shape Functions
    // --------------------------------------------------------------------------------------------

    //! @todo TT: fixed size dimensions could be used here, e.g.
    //! Eigen::Matrix<double, TDim,     Eigen::Dynamic> if B represents a normal gradient or
    //! Eigen::Matrix<double, VoigtDim, Eigen::Dynamic> if B calculates the engineering strain.
    //! However, it is harder to put them into one container. Maybe shared_ptr<Eigen::MatrixBase>?
    std::map<Node::eDof, Eigen::MatrixXd> mB;
    std::map<Node::eDof, const Eigen::MatrixXd*> mN;

    // Misc
    // --------------------------------------------------------------------------------------------

    double mDetJxWeightIPxSection;
    double mDetJacobian = 0;
    double mTotalMass = 0;

};

} /* namespace NuTo */
