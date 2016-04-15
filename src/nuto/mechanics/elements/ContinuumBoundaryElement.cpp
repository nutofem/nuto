/*
 * ContinuumBoundaryElementBase.cpp
 *
 *  Created on: 5 Mar 2016
 *      Author: vhirtham
 */

#include "nuto/mechanics/elements/ContinuumBoundaryElement.h"
#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/EvaluateDataContinuumBoundary.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/structures/StructureBase.h"

template <int TDim>
NuTo::ContinuumBoundaryElement<TDim>::ContinuumBoundaryElement(const ContinuumElement<TDim> *rBaseElement, int rSurfaceId)
: ElementBase::ElementBase(rBaseElement->GetStructure(), rBaseElement->GetElementDataType(), rBaseElement->GetIpDataType(0), rBaseElement->GetInterpolationType()),
  mBaseElement(rBaseElement),
  mSurfaceId(rSurfaceId),
  mBoundaryConditionType(BoundaryType::NOT_SET)
{}

template <int TDim>
NuTo::Error::eError NuTo::ContinuumBoundaryElement<TDim>::Evaluate(const ConstitutiveInputMap& rInput, std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput)
{
    EvaluateDataContinuumBoundary<TDim> data;
    ExtractAllNecessaryDofValues(data);

    auto constitutiveOutput = GetConstitutiveOutputMap(rElementOutput, data);
    auto constitutiveInput  = GetConstitutiveInputMap(constitutiveOutput, data);
    constitutiveInput.insert(rInput.begin(), rInput.end());

    for (int theIP = 0; theIP < GetNumIntegrationPoints(); theIP++)
    {
        CalculateNMatrixBMatrixDetJacobian(data, theIP);
        CalculateConstitutiveInputs(constitutiveInput, data);

        try
        {
            ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(theIP);
            for(auto itOutput : constitutiveOutput)
                if(itOutput.second!=nullptr) //check nullptr because of static data
                    itOutput.second->SetIsCalculated(false);
            Error::eError error = constitutivePtr->Evaluate<TDim>(this, theIP, constitutiveInput, constitutiveOutput);
            if (error != Error::SUCCESSFUL)
                return error;            
            CalculateGradientDamageBoundaryConditionParameters(data, constitutivePtr);
            for(auto itOutput : constitutiveOutput)
                if(itOutput.second!=nullptr && !itOutput.second->GetIsCalculated()) //check nullptr because of static data
                    throw MechanicsException(__PRETTY_FUNCTION__,std::string("Output ")+Constitutive::OutputToString(itOutput.first)+" not calculated by constitutive law");

        } catch (NuTo::MechanicsException& e)
        {
            e.AddMessage(__PRETTY_FUNCTION__, "error evaluating the constitutive model.");
            throw e;
        }
        CalculateElementOutputs(rElementOutput, data, theIP);
    }
    return Error::SUCCESSFUL;
}

template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::ExtractAllNecessaryDofValues(EvaluateDataContinuumBoundary<TDim> &rData)
{
    // needs optimization,
    // not all dofs might be needed...

    const std::set<Node::eDof>& dofs = mInterpolationType->GetDofs();
    for (auto dof : dofs)
        if (mInterpolationType->IsConstitutiveInput(dof))
            rData.mNodalValues[dof] = mBaseElement->ExtractNodeValues(0, dof);

    rData.mNodalValues[Node::COORDINATES] = mBaseElement->ExtractNodeValues(0, Node::COORDINATES);

    if (mStructure->GetNumTimeDerivatives() >= 1)
        for (auto dof : dofs)
            if (mInterpolationType->IsConstitutiveInput(dof))
                rData.mNodalValues_dt1[dof] = mBaseElement->ExtractNodeValues(1, dof);
}


template<int TDim>
NuTo::ConstitutiveOutputMap NuTo::ContinuumBoundaryElement<TDim>::GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
                                                                                   EvaluateDataContinuumBoundary<TDim>& rData) const
{
    ConstitutiveOutputMap constitutiveOutput;

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::INTERNAL_GRADIENT:
            FillConstitutiveOutputMapInternalGradient(constitutiveOutput, it.second->GetBlockFullVectorDouble(), rData);
            break;

        case Element::HESSIAN_0_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian0(constitutiveOutput, it.second->GetBlockFullMatrixDouble(), rData);
            break;

        case Element::HESSIAN_1_TIME_DERIVATIVE:
            FillConstitutiveOutputMapHessian1(constitutiveOutput, it.second->GetBlockFullMatrixDouble(), rData);
            break;

        case Element::HESSIAN_2_TIME_DERIVATIVE:
            throw MechanicsException(__PRETTY_FUNCTION__, "Case not handled! IMPORTANT: Everything must be set to zero here!!!");
            break;

        case Element::UPDATE_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
            break;

        case Element::UPDATE_TMP_STATIC_DATA:
            constitutiveOutput[NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA] = 0;
            break;

        case Element::IP_DATA:
            FillConstitutiveOutputMapIpData(constitutiveOutput, it.second->GetIpData(), rData);
            break;

        case Element::GLOBAL_ROW_DOF:
            mBaseElement->CalculateGlobalRowDofs(it.second->GetBlockFullVectorInt());
            break;

        case Element::GLOBAL_COLUMN_DOF:
            mBaseElement->CalculateGlobalColumnDofs(it.second->GetBlockFullVectorInt());
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element  output not implemented.");
        }
    }
    return constitutiveOutput;
}


template<int TDim>
NuTo::ConstitutiveInputMap NuTo::ContinuumBoundaryElement<TDim>::GetConstitutiveInputMap(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                                                EvaluateDataContinuumBoundary<TDim>& rData) const
{
    ConstitutiveInputMap constitutiveInputMap =  GetConstitutiveLaw(0)->GetConstitutiveInputs(rConstitutiveOutput, *GetInterpolationType());

    for (auto& itInput : constitutiveInputMap)
    {
        switch (itInput.first)
        {
        case Constitutive::Input::ENGINEERING_STRAIN:
            itInput.second = &(rData.mEngineeringStrain);
            break;

        case Constitutive::Input::NONLOCAL_EQ_STRAIN:
            itInput.second = &(rData.mNonlocalEqStrain);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY:
            itInput.second = &(rData.mRelativeHumidity);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION:
            itInput.second = &(rData.mWaterVolumeFraction);
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive input " + Constitutive::InputToString(itInput.first) + " cannot be calculated by this element type.");
        }
    }
    return constitutiveInputMap;
}





template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateNMatrixBMatrixDetJacobian(EvaluateDataContinuumBoundary<TDim> &rData, int rTheIP) const
{

    const InterpolationBase& interpolationTypeCoords = mInterpolationType->Get(Node::COORDINATES);

    Eigen::Matrix<double,TDim-1,1>  ipCoordsSurface = CalculateIPCoordinatesSurface(rTheIP);
    Eigen::Matrix<double,TDim,1>    ipCoordsNatural = interpolationTypeCoords.CalculateNaturalSurfaceCoordinates(ipCoordsSurface, mSurfaceId);

    // #######################################
    // ##  Calculate the surface jacobian
    // ## = || [dX / dXi] * [dXi / dAlpha] ||
    // #######################################
    Eigen::MatrixXd derivativeShapeFunctionsNatural     = interpolationTypeCoords.CalculateDerivativeShapeFunctionsNatural(ipCoordsNatural);
    const Eigen::Matrix<double,TDim,TDim> jacobian      = mBaseElement->CalculateJacobian(derivativeShapeFunctionsNatural, rData.mNodalValues[Node::COORDINATES]);// = [dX / dXi]

    const Eigen::MatrixXd derivativeNaturalSurfaceCoordinates   = interpolationTypeCoords.CalculateDerivativeNaturalSurfaceCoordinates(ipCoordsSurface, mSurfaceId); // = [dXi / dAlpha]
    rData.mDetJacobian = (jacobian * derivativeNaturalSurfaceCoordinates).norm();

    if (rData.mDetJacobian == 0)
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Determinant of the Jacobian is zero, no inversion possible.");
    }

    const Eigen::Matrix<double,TDim,TDim> invJacobian   = jacobian.inverse();


    for (auto dof : mInterpolationType->GetDofs())
    {
        if (dof == Node::COORDINATES)
            continue;
        const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
        rData.mN[dof] = interpolationType.CalculateMatrixN(ipCoordsNatural);

        rData.mB[dof] = mBaseElement->CalculateMatrixB(dof, interpolationType.CalculateDerivativeShapeFunctionsNatural(ipCoordsNatural), invJacobian);
    }
}



template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateConstitutiveInputs(const ConstitutiveInputMap& rConstitutiveInput, EvaluateDataContinuumBoundary<TDim> &rData)
{
    for (auto it : rConstitutiveInput)
    {
        switch (it.first)
        {
        case Constitutive::Input::ENGINEERING_STRAIN:
            rData.mEngineeringStrain.AsVector() = rData.mB.at(Node::DISPLACEMENTS) * rData.mNodalValues.at(Node::DISPLACEMENTS);
            break;

        case Constitutive::Input::NONLOCAL_EQ_STRAIN:
            rData.mNonlocalEqStrain.AsScalar() = rData.mN.at(Node::NONLOCALEQSTRAIN) * rData.mNodalValues.at(Node::NONLOCALEQSTRAIN);
            break;

        case Constitutive::Input::RELATIVE_HUMIDITY:
            rData.mRelativeHumidity.AsScalar() = rData.mN.at(Node::RELATIVEHUMIDITY) * rData.mNodalValues.at(Node::RELATIVEHUMIDITY);
            break;

        case Constitutive::Input::WATER_VOLUME_FRACTION:
            rData.mWaterVolumeFraction.AsScalar() = rData.mN.at(Node::WATERVOLUMEFRACTION) * rData.mNodalValues.at(Node::WATERVOLUMEFRACTION);
            break;

        case Constitutive::Input::TIME_STEP:
        case Constitutive::Input::CALCULATE_STATIC_DATA:
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive input for " + Constitutive::InputToString(it.first) + " not implemented.");
        }
    }
}

template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
                                                          EvaluateDataContinuumBoundary<TDim> &rData,
                                                          int rTheIP) const
{
    rData.mDetJxWeightIPxSection = CalculateDetJxWeightIPxSection(rData.mDetJacobian, rTheIP); // formerly known as "factor"

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::INTERNAL_GRADIENT:
            CalculateElementOutputInternalGradient(it.second->GetBlockFullVectorDouble(), rData, rTheIP);
            break;

        case Element::HESSIAN_0_TIME_DERIVATIVE:
            CalculateElementOutputHessian0(it.second->GetBlockFullMatrixDouble(), rData, rTheIP);
            break;

        case Element::HESSIAN_1_TIME_DERIVATIVE:
            break;

        case Element::HESSIAN_2_TIME_DERIVATIVE:
            break;

        case Element::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
            break;

        case Element::UPDATE_STATIC_DATA:
        case Element::UPDATE_TMP_STATIC_DATA:
            break;
        case Element::IP_DATA:
            CalculateElementOutputIpData(it.second->GetIpData(), rData, rTheIP);
            break;
        case Element::GLOBAL_ROW_DOF:
        case Element::GLOBAL_COLUMN_DOF:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
        }
    }
}


template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputInternalGradient(BlockFullVector<double>& rInternalGradient,
                                                                         EvaluateDataContinuumBoundary<TDim> &rData,
                                                                         int rTheIP) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::DISPLACEMENTS:
            break;

        case Node::NONLOCALEQSTRAIN:
        {
            if (rData.mBCType == BoundaryType::ROBIN_INHOMOGENEOUS)
                rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection / rData.mAlpha * rData.mN.at(Node::NONLOCALEQSTRAIN).transpose() * (rData.mNonlocalEqStrain[0] - rData.mLocalEqStrain[0]);

            break;
        }

        case Node::RELATIVEHUMIDITY:
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * rData.mInternalGradientRH_Boundary_N;
            break;

        case Node::WATERVOLUMEFRACTION:
            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * rData.mInternalGradientWV_Boundary_N;
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Element output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");
        }
    }
}



template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0,
                                                                 EvaluateDataContinuumBoundary<TDim> &rData,
                                                                 int rTheIP) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            auto& hessian0 = rHessian0(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::DISPLACEMENTS):
            case Node::CombineDofs(Node::DISPLACEMENTS, Node::NONLOCALEQSTRAIN):
                break;

            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::DISPLACEMENTS):
                if (rData.mBCType == BoundaryType::ROBIN_INHOMOGENEOUS)
                    hessian0 -= rData.mDetJxWeightIPxSection / rData.mAlpha * rData.mN.at(dofRow).transpose() * rData.mTangentLocalEqStrainStrain.transpose() * rData.mB.at(dofCol);
                break;

            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::NONLOCALEQSTRAIN):
                if (rData.mBCType == BoundaryType::ROBIN_INHOMOGENEOUS)
                    hessian0 += rData.mN.at(dofRow).transpose() * rData.mN.at(dofRow) * rData.mDetJxWeightIPxSection / rData.mAlpha;
                break;

            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
                hessian0 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * rData.mInternalGradientRH_dRH_Boundary_NN_H0 * rData.mN.at(dofCol);
                break;

            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
                hessian0 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * rData.mInternalGradientWV_dWV_Boundary_NN_H0 * rData.mN.at(dofCol);
                break;

            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
                break;
            default:
                /*******************************************************\
                |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
                \*******************************************************/
                 case Node::CombineDofs(Node::DISPLACEMENTS,         Node::RELATIVEHUMIDITY):
                 case Node::CombineDofs(Node::DISPLACEMENTS,         Node::WATERVOLUMEFRACTION):
                 case Node::CombineDofs(Node::RELATIVEHUMIDITY,      Node::DISPLACEMENTS):
                 case Node::CombineDofs(Node::WATERVOLUMEFRACTION,   Node::DISPLACEMENTS):
                     continue;
                throw MechanicsException(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}

template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateElementOutputIpData(ElementOutputIpData& rIpData, EvaluateDataContinuumBoundary<TDim> &rData, int rTheIP) const
{
    for (auto& it : rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    {                                       // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::ENGINEERING_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mEngineeringStrainVisualize);
            break;
        case NuTo::IpData::ENGINEERING_STRESS:
            it.second.col(rTheIP) = std::move(rData.mEngineeringStressVisualize);
            break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mEngineeringPlasticStrainVisualize);
            break;
        case NuTo::IpData::DAMAGE:
            it.second.col(rTheIP) = std::move(rData.mDamage);
            break;
        case NuTo::IpData::EXTRAPOLATION_ERROR:
            it.second.col(rTheIP) = std::move(rData.mExtrapolationError);
            break;
        case NuTo::IpData::LOCAL_EQ_STRAIN:
            it.second.col(rTheIP) = std::move(rData.mLocalEqStrain);
            break;
        default:
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Ip data not implemented.");
        }
    }
}

template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapInternalGradient(ConstitutiveOutputMap& rConstitutiveOutput,
                                                                            BlockFullVector<double>& rInternalGradient,
                                                                            EvaluateDataContinuumBoundary<TDim>& rData) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        rInternalGradient[dofRow].Resize(mInterpolationType->Get(dofRow).GetNumDofs());
        switch (dofRow)
        {
        case Node::DISPLACEMENTS:
            break;

        case Node::NONLOCALEQSTRAIN:
            rConstitutiveOutput[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &(rData.mLocalEqStrain);
            rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &(rData.mNonlocalParameterXi);
            break;

        case Node::RELATIVEHUMIDITY:
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N] = &(rData.mInternalGradientRH_Boundary_N);
            break;

        case Node::WATERVOLUMEFRACTION:
            rConstitutiveOutput[NuTo::Constitutive::Output::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N] = &(rData.mInternalGradientWV_Boundary_N);
            break;

        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output INTERNAL_GRADIENT for " + Node::DofToString(dofRow) + " not implemented.");

        }
    }
}

template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapHessian0(ConstitutiveOutputMap& rConstitutiveOutput,
                                                                    BlockFullMatrix<double>& rHessian0,
                                                                    EvaluateDataContinuumBoundary<TDim> &rData) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& dofSubMatrix = rHessian0(dofRow, dofCol);
            dofSubMatrix.Resize(mInterpolationType->Get(dofRow).GetNumDofs(), mInterpolationType->Get(dofCol).GetNumDofs());
            dofSubMatrix.setZero();

            switch (Node::CombineDofs(dofRow, dofCol))
            {

            case Node::CombineDofs(Node::DISPLACEMENTS, Node::DISPLACEMENTS):
            case Node::CombineDofs(Node::DISPLACEMENTS, Node::NONLOCALEQSTRAIN):
                break;

            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::DISPLACEMENTS):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN] = &(rData.mTangentLocalEqStrainStrain);
                break;

            case Node::CombineDofs(Node::NONLOCALEQSTRAIN, Node::NONLOCALEQSTRAIN):
                rConstitutiveOutput[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &(rData.mNonlocalParameterXi);
                break;

            case Node::CombineDofs(Node::RELATIVEHUMIDITY, Node::RELATIVEHUMIDITY):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0] = &rData.mInternalGradientRH_dRH_Boundary_NN_H0;
                break;

            case Node::CombineDofs(Node::WATERVOLUMEFRACTION, Node::WATERVOLUMEFRACTION):
                rConstitutiveOutput[NuTo::Constitutive::Output::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0] = &rData.mInternalGradientWV_dWV_Boundary_NN_H0;
                break;


           /*******************************************************\
           |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
           \*******************************************************/
            case Node::CombineDofs(Node::DISPLACEMENTS, Node::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::DISPLACEMENTS, Node::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::RELATIVEHUMIDITY, Node::DISPLACEMENTS):
            case Node::CombineDofs(Node::RELATIVEHUMIDITY, Node::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::WATERVOLUMEFRACTION, Node::DISPLACEMENTS):
            case Node::CombineDofs(Node::WATERVOLUMEFRACTION, Node::RELATIVEHUMIDITY):
            {
                continue;
            }
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_0_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}



template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap& rConstitutiveOutput,
                                                                    BlockFullMatrix<double>& rHessian0,
                                                                    EvaluateDataContinuumBoundary<TDim> &rData) const
{
    for (auto dofRow : mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : mInterpolationType->GetActiveDofs())
        {
            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& dofSubMatrix = rHessian0(dofRow, dofCol);
            dofSubMatrix.Resize(mInterpolationType->Get(dofRow).GetNumDofs(), mInterpolationType->Get(dofCol).GetNumDofs());
            dofSubMatrix.setZero();

            switch (Node::CombineDofs(dofRow, dofCol))
            {



           /*******************************************************\
           |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
           \*******************************************************/
            case Node::CombineDofs(Node::DISPLACEMENTS,         Node::DISPLACEMENTS):
            case Node::CombineDofs(Node::DISPLACEMENTS,         Node::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::DISPLACEMENTS,         Node::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::RELATIVEHUMIDITY,      Node::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::WATERVOLUMEFRACTION,   Node::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::RELATIVEHUMIDITY,      Node::DISPLACEMENTS):
            case Node::CombineDofs(Node::RELATIVEHUMIDITY,      Node::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::WATERVOLUMEFRACTION,   Node::DISPLACEMENTS):
            case Node::CombineDofs(Node::WATERVOLUMEFRACTION,   Node::RELATIVEHUMIDITY):
                continue;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive output HESSIAN_1_TIME_DERIVATIVE for "
                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
            }
        }
    }
}


//template<int TDim>
//void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian2, EvaluateDataContinuum<TDim> &rData) const
//{
//    for (auto dofRow : mInterpolationType->GetActiveDofs())
//    {
//        for (auto dofCol : mInterpolationType->GetActiveDofs())
//        {
//            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& dofSubMatrix = rHessian2(dofRow, dofCol);
//            dofSubMatrix.Resize(mInterpolationType->Get(dofRow).GetNumDofs(), mInterpolationType->Get(dofCol).GetNumDofs());
//            dofSubMatrix.setZero();

//            switch (Node::CombineDofs(dofRow, dofCol))
//            {
//            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
//                break;
//            default:
//                throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] Constitutive output HESSIAN_2_TIME_DERIVATIVE for "
//                        "(" + Node::DofToString(dofRow) + "," + Node::DofToString(dofCol) + ") not implemented.");
//            }
//        }
//    }
//}
//
template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::FillConstitutiveOutputMapIpData(ConstitutiveOutputMap& rConstitutiveOutput, ElementOutputIpData& rIpData, EvaluateDataContinuumBoundary<TDim> &rData) const
{

    for (auto& it : rIpData.GetIpDataMap()) // this reference here is _EXTREMLY_ important, since the GetIpDataMap() contains a
    {                                       // FullMatrix VALUE and you want to access this value by reference. Without the &, a tmp copy would be made.
        switch (it.first)
        {
        case NuTo::IpData::ENGINEERING_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE] = &(rData.mEngineeringStrainVisualize);
            break;
        case NuTo::IpData::ENGINEERING_STRESS:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE] = &(rData.mEngineeringStressVisualize);
            break;
        case NuTo::IpData::ENGINEERING_PLASTIC_STRAIN:
            it.second.Resize(6, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::ENGINEERING_PLASTIC_STRAIN_VISUALIZE] = &(rData.mEngineeringPlasticStrainVisualize);
            break;
        case NuTo::IpData::DAMAGE:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::DAMAGE] = &(rData.mDamage);
            break;
        case NuTo::IpData::LOCAL_EQ_STRAIN:
            it.second.Resize(1, GetNumIntegrationPoints());
            rConstitutiveOutput[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &(rData.mLocalEqStrain);
            break;
        default:
            throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "] this ip data type is not implemented.");
        }
    }
}





template <int TDim>
const Eigen::Vector3d NuTo::ContinuumBoundaryElement<TDim>::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    Eigen::VectorXd naturalSurfaceIpCoordinates;
    switch (GetStructure()->GetDimension())
    {
        case 1:
        {
            double ipCoordinate;
            GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, ipCoordinate);
            naturalSurfaceIpCoordinates.resize(1);
            naturalSurfaceIpCoordinates(0) = ipCoordinate;
            break;
        }
        case 2:
        {
            double ipCoordinates[2];
            GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, ipCoordinates);
            naturalSurfaceIpCoordinates.resize(2);
            naturalSurfaceIpCoordinates(0) = ipCoordinates[0];
            naturalSurfaceIpCoordinates(1) = ipCoordinates[1];
            break;
        }
        case 3:
        {
            double ipCoordinates[3];
            GetIntegrationType()->GetLocalIntegrationPointCoordinates3D(rIpNum, ipCoordinates);
            naturalSurfaceIpCoordinates.resize(3);
            naturalSurfaceIpCoordinates(0) = ipCoordinates[0];
            naturalSurfaceIpCoordinates(1) = ipCoordinates[1];
            naturalSurfaceIpCoordinates(2) = ipCoordinates[2];
            break;
        }
        default:
            break;
    }

    Eigen::VectorXd naturalIpCoordinates = mInterpolationType->Get(Node::COORDINATES).CalculateNaturalSurfaceCoordinates(naturalSurfaceIpCoordinates, mSurfaceId);

    Eigen::VectorXd matrixN = mInterpolationType->Get(Node::COORDINATES).CalculateMatrixN(naturalIpCoordinates);
    Eigen::VectorXd nodeCoordinates = ExtractNodeValues(0, Node::COORDINATES);

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();
    globalIntegrationPointCoordinates.segment(0, GetLocalDimension()) = matrixN * nodeCoordinates;

    return globalIntegrationPointCoordinates;
}


#ifdef ENABLE_VISUALIZE
template <int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::Visualize(VisualizeUnstructuredGrid& rVisualize, const std::list<std::shared_ptr<NuTo::VisualizeComponent>>& rVisualizationList)
{
    if (GetStructure()->GetVerboseLevel() > 10)
        std::cout << __PRETTY_FUNCTION__ << "Pleeeaaase, implement the visualization for me!!!" << std::endl;
}
#endif


template<int TDim>
void NuTo::ContinuumBoundaryElement<TDim>::CalculateGradientDamageBoundaryConditionParameters(EvaluateDataContinuumBoundary<TDim>& rData, const ConstitutiveBase* rConstitutiveLaw) const
{
    //VHIRTHAMTODO Temporary solution --- find better one ---> Problem: If dof NONLOCALEQSTRAIN does not exists, we will get an exception.
    const std::set<Node::eDof>& DofTypes = mStructure->GetDofStatus().GetDofTypes();
    if(DofTypes.find(Node::NONLOCALEQSTRAIN)==DofTypes.end())
        return;

    if (not mInterpolationType->IsActive(Node::NONLOCALEQSTRAIN))
        return;

    rData.mAlpha = std::sqrt(rData.mNonlocalParameterXi[0]);

    double e0 = rConstitutiveLaw->GetParameterDouble(Constitutive::eConstitutiveParameter::TENSILE_STRENGTH) /
                rConstitutiveLaw->GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS);
    if (mBoundaryConditionType == BoundaryType::MACAULAY)
    {
        // determine state ...
        bool switchToNeumann = (rData.mLocalEqStrain[0] > rData.mNonlocalEqStrain[0]) and (rData.mLocalEqStrain[0] > e0);

        // ... and use existing implementations
        if (switchToNeumann)
        {
            rData.mBCType = BoundaryType::NEUMANN_HOMOGENEOUS;
//            std::cout << "Macaulay Culcin helps out in element " << GetStructure()->ElementGetId(this) << std::endl;

        }
        else
            rData.mBCType = BoundaryType::ROBIN_INHOMOGENEOUS;
    }
    else
    {
        rData.mBCType = mBoundaryConditionType;
    }
}


namespace NuTo
{
template<>
NuTo::ConstitutiveStaticDataBase* ContinuumBoundaryElement<1>::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticData1D(this);
}
template<>
NuTo::ConstitutiveStaticDataBase* ContinuumBoundaryElement<2>::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticData2D(this);
}
template<>
NuTo::ConstitutiveStaticDataBase* ContinuumBoundaryElement<3>::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticData3D(this);
}






template<>
Eigen::Matrix<double,0,1>  ContinuumBoundaryElement<1>::CalculateIPCoordinatesSurface(int rTheIP) const
{
    return Eigen::Matrix<double,0,1>();
}

template<>
Eigen::Matrix<double,1,1>  ContinuumBoundaryElement<2>::CalculateIPCoordinatesSurface(int rTheIP) const
{
    double tmp;
    GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rTheIP, tmp);
    Eigen::Matrix<double,1,1> ipCoordinatesSurface;
    ipCoordinatesSurface(0) = tmp;
    return ipCoordinatesSurface;
}

template<>
Eigen::Matrix<double,2,1>  ContinuumBoundaryElement<3>::CalculateIPCoordinatesSurface(int rTheIP) const
{
    double tmp[2];
    GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rTheIP, tmp);
    Eigen::Matrix<double,2,1> ipCoordinatesSurface;
    ipCoordinatesSurface(0) = tmp[0];
    ipCoordinatesSurface(1) = tmp[1];
    return ipCoordinatesSurface;
}



template<>
double NuTo::ContinuumBoundaryElement<1>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{

    return  mBaseElement->mSection->GetArea();
}

template<>
double NuTo::ContinuumBoundaryElement<2>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return  rDetJacobian *
            mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP) *
            mBaseElement->mSection->GetThickness();

}

template<>
double NuTo::ContinuumBoundaryElement<3>::CalculateDetJxWeightIPxSection(double rDetJacobian, int rTheIP) const
{
    return  rDetJacobian *
            mElementData->GetIntegrationType()->GetIntegrationPointWeight(rTheIP);
}


template<>
const ContinuumBoundaryElement<1>& ContinuumBoundaryElement<1>::AsContinuumBoundaryElement1D() const
{
    return *this;
}

template<>
const ContinuumBoundaryElement<2>& ContinuumBoundaryElement<2>::AsContinuumBoundaryElement2D() const
{
    return *this;
}

template<>
const ContinuumBoundaryElement<3>& ContinuumBoundaryElement<3>::AsContinuumBoundaryElement3D() const
{
    return *this;
}

template<>
ContinuumBoundaryElement<1>& ContinuumBoundaryElement<1>::AsContinuumBoundaryElement1D()
{
    return *this;
}

template<>
ContinuumBoundaryElement<2>& ContinuumBoundaryElement<2>::AsContinuumBoundaryElement2D()
{
    return *this;
}

template<>
ContinuumBoundaryElement<3>& ContinuumBoundaryElement<3>::AsContinuumBoundaryElement3D()
{
    return *this;
}

} //namespace NuTo

template class NuTo::ContinuumBoundaryElement<1>;
template class NuTo::ContinuumBoundaryElement<2>;
template class NuTo::ContinuumBoundaryElement<3>;
