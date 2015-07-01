/*
 * BoundaryElement1D.cpp
 *
 *  Created on: 4 Jun 2015
 *      Author: ttitsche
 */

#include "nuto/mechanics/elements/BoundaryElement1D.h"
#include "nuto/mechanics/elements/Element1D.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/sections/SectionBase.h"

#include "nuto/mechanics/constitutive/ConstitutiveTangentLocal.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/constitutive/moistureTransport/RelativeHumidity.h"
#include "nuto/mechanics/constitutive/moistureTransport/WaterVolumeFraction.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/DeformationGradient1D.h"

#include "nuto/mechanics/constitutive/mechanics/NonlocalEqStrain.h"
#include "nuto/mechanics/constitutive/mechanics/LocalEqStrain.h"

NuTo::BoundaryElement1D::BoundaryElement1D(const ElementBase* rBaseElement, int rSurfaceId) :
    NuTo::BoundaryElementBase::BoundaryElementBase(rBaseElement, rSurfaceId), mBoundaryRelativeHumidity(0.0), mBoundaryWaterVolumeFraction(0.0)
{

}

NuTo::Error::eError NuTo::BoundaryElement1D::Evaluate(boost::ptr_multimap<NuTo::Element::eOutput, NuTo::ElementOutputBase>& rElementOutput)
{
    try
    {
        const SectionBase* section(GetSection());
        if (section == 0)
            throw MechanicsException("[NuTo::BoundaryElement1D::Evaluate] no section allocated for element.");

        const std::set<Node::eAttributes>& dofs = mInterpolationType->GetDofs();
        const std::set<Node::eAttributes>& activeDofs = mInterpolationType->GetActiveDofs();

        int numActiveDofs = mInterpolationType->GetNumActiveDofs();

        // extract all node values and store them
        std::map<Node::eAttributes, Eigen::MatrixXd> nodalValues;
        for (auto dof : dofs)
        {
            nodalValues[dof] = ExtractNodeValues(0, dof);
        }

        // Gradient damage model
        DeformationGradient1D deformationGradient;
        NonlocalEqStrain nonlocalEqStrain;
        LocalEqStrain localEqStrain;
        ConstitutiveTangentLocal<1, 1> tangentLocalEqStrainStrain;
        ConstitutiveTangentLocal<1, 1> nonlocalParameter;

        //allocate relative humidity
        RelativeHumidity relativeHumidity;
        RelativeHumidity relativeHumidityD1;
        RelativeHumidity relativeHumidityBoundary;

        //allocate water phase fraction
        WaterVolumeFraction waterVolumeFraction;
        WaterVolumeFraction waterVolumeFractionD1;
        WaterVolumeFraction waterVolumeFractionBoundary;


        ConstitutiveTangentLocal<1,1> residualBoundarySurfaceWaterPhase;
        ConstitutiveTangentLocal<1,1> residualBoundarySurfaceVaporPhase;
        ConstitutiveTangentLocal<1,1> tangentSurfaceRelativeHumidityTransportCoefficient;
        ConstitutiveTangentLocal<1,1> tangentSurfaceWaterVolumeFractionTransportCoefficient;

        //allocate input list and output list
        std::map<NuTo::Constitutive::Input::eInput, const ConstitutiveInputBase*> constitutiveInputList;
        std::map<NuTo::Constitutive::Output::eOutput, ConstitutiveOutputBase*> constitutiveOutputList;


        //define constitutive input list
        for (auto dof : dofs)
        {
            if (mInterpolationType->IsConstitutiveInput(dof) == false)
                continue;
            switch (dof)
            {
            case Node::NONLOCALEQSTRAIN:
                constitutiveInputList[NuTo::Constitutive::Input::DEFORMATION_GRADIENT_1D] = &deformationGradient;
                constitutiveInputList[NuTo::Constitutive::Input::NONLOCAL_EQ_STRAIN] = &(nonlocalEqStrain);
                break;
            case Node::RELATIVEHUMIDITY:
            {
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY]             = &relativeHumidity;
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_BOUNDARY]    = &relativeHumidityBoundary;
                constitutiveInputList[NuTo::Constitutive::Input::RELATIVE_HUMIDITY_D1]          = &relativeHumidityD1;
                break;
            }
            case Node::WATERVOLUMEFRACTION:
            {
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION]         = &waterVolumeFraction;
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_BOUNDARY]= &waterVolumeFractionBoundary;
                constitutiveInputList[NuTo::Constitutive::Input::WATER_VOLUME_FRACTION_D1]      = &waterVolumeFractionD1;
                break;
            }
            default:
                break;
            }
        }



        //define constitutive output list
        for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
        {
            switch (it->first)
            {
            case Element::INTERNAL_GRADIENT:
                it->second->GetFullVectorDouble().Resize(numActiveDofs);
                //if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
                //on the global level
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::NONLOCALEQSTRAIN:
                        constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &localEqStrain;
                        constitutiveOutputList[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &nonlocalParameter;
                        break;

                    case Node::RELATIVEHUMIDITY:
                    {
                        if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::BOUNDARY_SURFACE_VAPOR_PHASE_RESIDUAL]  = &residualBoundarySurfaceVaporPhase;
                        }
                        break;
                    }
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::BOUNDARY_SURFACE_WATER_PHASE_RESIDUAL]  = &residualBoundarySurfaceWaterPhase;
                        }
                        break;
                    }

                    default:
                        break;
                    }
                }
                break;
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->GetFullMatrixDouble().setZero();
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                for (auto dof : activeDofs)
                {
                    switch (dof)
                    {
                    case Node::NONLOCALEQSTRAIN:
                        constitutiveOutputList[NuTo::Constitutive::Output::LOCAL_EQ_STRAIN] = &localEqStrain;
                        constitutiveOutputList[NuTo::Constitutive::Output::D_LOCAL_EQ_STRAIN_D_STRAIN_1D] = &tangentLocalEqStrainStrain;
                        constitutiveOutputList[NuTo::Constitutive::Output::NONLOCAL_PARAMETER_XI] = &nonlocalParameter;
                        break;
                    case Node::RELATIVEHUMIDITY:
                    {
                        if (activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::BOUNDARY_SURFACE_RELATIVE_HUMIDIY_TRANSPORT_COEFFICIENT]  = &tangentSurfaceRelativeHumidityTransportCoefficient;
                        }
                        break;
                    }
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if (activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            constitutiveOutputList[NuTo::Constitutive::Output::BOUNDARY_SURFACE_WATER_VOLUME_FRACTION_TRANSPORT_COEFFICIENT]  = &tangentSurfaceWaterVolumeFractionTransportCoefficient;
                        }
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
                break;
            case Element::HESSIAN_1_TIME_DERIVATIVE:
                it->second->GetFullMatrixDouble().Resize(numActiveDofs, numActiveDofs);
                it->second->SetSymmetry(true);
                it->second->SetConstant(true);
                break;
            case Element::UPDATE_STATIC_DATA:
                constitutiveOutputList[NuTo::Constitutive::Output::UPDATE_STATIC_DATA] = 0;
                break;
            case Element::GLOBAL_ROW_DOF:
            {
                const Eigen::VectorXi& globalRowDofsEigen = mBaseElement->AsElement1D()->CalculateGlobalRowDofs();
                std::vector<int> globalRowDofsStd(globalRowDofsEigen.data(), globalRowDofsEigen.data() + globalRowDofsEigen.rows());
                it->second->GetVectorInt() = globalRowDofsStd;
            }
                break;
            case Element::GLOBAL_COLUMN_DOF:
            {
                const Eigen::VectorXi& globalColumnDofsEigen = mBaseElement->AsElement1D()->CalculateGlobalColumnDofs();
                std::vector<int> globalColumnDofsStd(globalColumnDofsEigen.data(), globalColumnDofsEigen.data() + globalColumnDofsEigen.rows());
                it->second->GetVectorInt() = globalColumnDofsStd;
            }
                break;
            default:
                break;
            }
        }               //end for: constitutive output list

        std::map<Node::eAttributes, Eigen::VectorXd> shapeFunctions;
        std::map<Node::eAttributes, Eigen::MatrixXd> derivativeShapeFunctions;

        double factor = section->GetArea();

        Eigen::VectorXd dummy; // not needed in 1D
        Eigen::VectorXd naturalIpCoordinatesSurface = mInterpolationType->Get(Node::COORDINATES).CalculateNaturalSurfaceCoordinates(dummy, mSurfaceId);

        //        std::cout << "SurfaceID " << mSurfaceId << std::endl;
        //        std::cout << "naturalIpCoordinatesSurface " << naturalIpCoordinatesSurface.transpose() << std::endl;

        const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural = mInterpolationType->Get(Node::COORDINATES).CalculateDerivativeShapeFunctionsNatural(naturalIpCoordinatesSurface);
        double detJacobian = mBaseElement->AsElement1D()->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, nodalValues[Node::COORDINATES]);

        // calculate shape functions and their derivatives
        for (auto dof : dofs)
        {
            if (dof == Node::COORDINATES)
                continue;
            const InterpolationBase& interpolationType = mInterpolationType->Get(dof);
            shapeFunctions[dof] = interpolationType.CalculateShapeFunctions(naturalIpCoordinatesSurface);
            derivativeShapeFunctions[dof] = interpolationType.CalculateDerivativeShapeFunctionsNatural(naturalIpCoordinatesSurface) / detJacobian;
        }

        // define constitutive inputs
        for (auto dof : dofs)
        {
            if (mInterpolationType->IsConstitutiveInput(dof) == false)
                continue;
            switch (dof)
            {
            case Node::NONLOCALEQSTRAIN:
            {
                deformationGradient = mBaseElement->AsElement1D()->CalculateDeformationGradient(derivativeShapeFunctions.at(Node::DISPLACEMENTS), nodalValues.at(Node::DISPLACEMENTS));
                nonlocalEqStrain(0, 0) = (nodalValues[Node::NONLOCALEQSTRAIN] * shapeFunctions[Node::NONLOCALEQSTRAIN])(0, 0);
            }
                break;
            case Node::RELATIVEHUMIDITY:
            {
                relativeHumidity(0,0)               = (nodalValues[Node::RELATIVEHUMIDITY] * shapeFunctions[Node::RELATIVEHUMIDITY])(0,0);
                relativeHumidityD1(0,0)             = (ExtractNodeValues(1, Node::RELATIVEHUMIDITY) * shapeFunctions[Node::RELATIVEHUMIDITY])(0,0);
                relativeHumidityBoundary(0,0)       = mBoundaryRelativeHumidity;
            }
                break;
            case Node::WATERVOLUMEFRACTION:
            {
                waterVolumeFraction(0,0)            = (nodalValues[Node::WATERVOLUMEFRACTION] * shapeFunctions[Node::WATERVOLUMEFRACTION])(0,0);
                waterVolumeFractionD1(0,0)          = (ExtractNodeValues(1, Node::WATERVOLUMEFRACTION) * shapeFunctions[Node::WATERVOLUMEFRACTION])(0,0);
                waterVolumeFractionBoundary(0,0)    = mBoundaryWaterVolumeFraction;
            }
                break;
            default:
                break;
            }
        }

        ConstitutiveBase* constitutivePtr = GetConstitutiveLaw(0);
        try
        {
            Error::eError error = constitutivePtr->Evaluate1D(this, 0, constitutiveInputList, constitutiveOutputList);
            if (error != Error::SUCCESSFUL)
                return error;
        } catch (NuTo::MechanicsException &e)
        {
            e.AddMessage("[NuTo::BoundaryElement1D::Evaluate] error evaluating the constitutive model.");
            throw e;
        }

        //        std::cout << "Nonlocal parameter " << nonlocalParameter.GetValue(0) << std::endl;

        int numDofsNonlocalEqStrain = 0;
        NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> KkkMod;
        BoundaryType::eType currentBCType = mBoundaryConditionType;
        double alpha = 0.;

        if (activeDofs.find(Node::NONLOCALEQSTRAIN) != activeDofs.end())
        {
            numDofsNonlocalEqStrain = mInterpolationType->Get(Node::NONLOCALEQSTRAIN).GetNumDofs();
            KkkMod.Resize(numDofsNonlocalEqStrain, numDofsNonlocalEqStrain);

            alpha = std::sqrt(nonlocalParameter.GetValue(0));

            //            std::cout << "local eq strain    " << localEqStrain[0] << std::endl;
            //            std::cout << "nonlocal eq strain " << nonlocalEqStrain[0] << std::endl;


            double e0 = constitutivePtr->GetTensileStrength() / constitutivePtr->GetYoungsModulus();
            if (mBoundaryConditionType == BoundaryType::MACAULAY)
            {


                // determine state ...
                bool switchToNeumann = (localEqStrain[0] > nonlocalEqStrain[0]) and (localEqStrain[0] > e0);

                // ... and use existing implementations
                if (switchToNeumann)
                {
                    currentBCType = BoundaryType::NEUMANN_HOMOGENEOUS;
                    std::cout << "Macaulay Culcin helps out in element " << mBaseElement->GetStructure()->ElementGetId(this) << std::endl;

                }
                else
                    currentBCType = BoundaryType::ROBIN_INHOMOGENEOUS;
            }

            switch (currentBCType)
            {
            case BoundaryType::NOT_SET:
                throw MechanicsException("[NuTo::BoundaryElement1D::Evaluate] Boundary condition type not set! ");
                break;
            case BoundaryType::ROBIN_INHOMOGENEOUS:
                KkkMod = shapeFunctions[Node::NONLOCALEQSTRAIN] * shapeFunctions[Node::NONLOCALEQSTRAIN].transpose() * factor / alpha;
                break;
            default:
                break;
            }

        }

        //calculate output
        for (auto it = rElementOutput.begin(); it != rElementOutput.end(); it++)
        {
            switch (it->first)
            {
            case Element::INTERNAL_GRADIENT:
            {
                //if the stiffness matrix is constant, the corresponding internal force is calculated via the Kd
                //on the global level
                for (auto dof : activeDofs)
                {
                    int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                    switch (dof)
                    {
                    case Node::NONLOCALEQSTRAIN:
                    {
                        // add Nt c eqStrain
                        NuTo::FullVector<double, Eigen::Dynamic>& f = it->second->GetFullVectorDouble();
                        if (currentBCType == BoundaryType::ROBIN_INHOMOGENEOUS)
                        {
                            f.block(startIndex, 0, numDofsNonlocalEqStrain, 1) += shapeFunctions[Node::NONLOCALEQSTRAIN] * factor / alpha * (nonlocalEqStrain[0] - localEqStrain[0]);
                        }
                        break;
                    }
                    case Node::RELATIVEHUMIDITY:
                    {
                        if(activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            mBaseElement-> AsElement1D() -> AddDetJNtX( shapeFunctions.at(Node::RELATIVEHUMIDITY),
                                                                        residualBoundarySurfaceVaporPhase,
                                                                        factor,
                                                                        startIndex,
                                                                        it->second->GetFullVectorDouble());
                        }
                        break;
                    }
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if(activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            mBaseElement-> AsElement1D() -> AddDetJNtX( shapeFunctions.at(Node::WATERVOLUMEFRACTION),
                                                                        residualBoundarySurfaceWaterPhase,
                                                                        factor,
                                                                        startIndex,
                                                                        it->second->GetFullVectorDouble());
                        }
                        break;
                    }
                    default:
                        break;
                    }
                }
            }
                break;
            case Element::HESSIAN_0_TIME_DERIVATIVE:
            {
                for (auto dof : activeDofs)
                {
                    int startIndex = mInterpolationType->Get(dof).GetLocalStartIndex();
                    switch (dof)
                    {
                    case Node::NONLOCALEQSTRAIN:
                        // add the modified Kkk to the element output
                        it->second->GetFullMatrixDouble().AddBlock(startIndex, startIndex, KkkMod);
                        if(currentBCType == BoundaryType::ROBIN_INHOMOGENEOUS)
                        {
                            Eigen::MatrixXd Ked = -tangentLocalEqStrainStrain(0)*factor*shapeFunctions[Node::NONLOCALEQSTRAIN]*derivativeShapeFunctions[Node::DISPLACEMENTS].transpose();
                            int dStart = mInterpolationType->Get(Node::DISPLACEMENTS).GetLocalStartIndex();
                            int eStart = startIndex;

                            it->second->GetFullMatrixDouble().AddBlock(eStart, dStart, Ked);
                        }
                        break;
                    case Node::RELATIVEHUMIDITY:
                    {
                        if(activeDofs.find(Node::WATERVOLUMEFRACTION) != activeDofs.end())
                        {
                            auto RelHumShapeFunction = shapeFunctions.at(Node::RELATIVEHUMIDITY);

                            mBaseElement-> AsElement1D() -> AddDetJNtXN(RelHumShapeFunction,
                                                                        RelHumShapeFunction,
                                                                        tangentSurfaceRelativeHumidityTransportCoefficient,
                                                                        factor,
                                                                        startIndex,
                                                                        startIndex,
                                                                        it->second->GetFullMatrixDouble());
                        }
                        break;
                    }
                    case Node::WATERVOLUMEFRACTION:
                    {
                        if(activeDofs.find(Node::RELATIVEHUMIDITY) != activeDofs.end())
                        {
                            auto WatVolShapeFunction = shapeFunctions.at(Node::WATERVOLUMEFRACTION);

                            mBaseElement-> AsElement1D() -> AddDetJNtXN(WatVolShapeFunction,
                                                                        WatVolShapeFunction,
                                                                        tangentSurfaceWaterVolumeFractionTransportCoefficient,
                                                                        factor,
                                                                        startIndex,
                                                                        startIndex,
                                                                        it->second->GetFullMatrixDouble());
                        }
                        break;
                    }
                    default:
                        break;
                    }

                }
            }
                break;
            default:
                break;
            }
        }

    } catch (NuTo::MechanicsException& e)
    {
        std::stringstream ss;
        ss << mBaseElement->GetStructure()->ElementGetId(this);
        e.AddMessage("[NuTo::BoundaryElement1D::Evaluate] Error evaluating element data of element " + ss.str() + ".");
        throw e;
    }
    return Error::SUCCESSFUL;

}

double NuTo::BoundaryElement1D::GetBoundaryRelativeHumidity() const
{
    return mBoundaryRelativeHumidity;
}

void NuTo::BoundaryElement1D::SetBoundaryRelativeHumidity(double rBoundaryRelativeHumidity)
{
    mBoundaryRelativeHumidity = rBoundaryRelativeHumidity;
}

double NuTo::BoundaryElement1D::GetBoundaryWaterVolumeFraction() const
{
    return mBoundaryWaterVolumeFraction;
}

NuTo::ConstitutiveStaticDataBase* NuTo::BoundaryElement1D::AllocateStaticData(const ConstitutiveBase* rConstitutiveLaw) const
{
    return rConstitutiveLaw->AllocateStaticDataEngineeringStress_EngineeringStrain1D(this);
}

void NuTo::BoundaryElement1D::SetBoundaryWaterVolumeFraction(double rBoundaryWaterVolumeFraction)
{
    mBoundaryWaterVolumeFraction = rBoundaryWaterVolumeFraction;
}

//! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
const NuTo::BoundaryElement1D* NuTo::BoundaryElement1D::AsBoundaryElement1D() const
{
    return this;
}

//! @brief cast the base pointer to an BoundaryElement1D, otherwise throws an exception
NuTo::BoundaryElement1D* NuTo::BoundaryElement1D::AsBoundaryElement1D()
{
    return this;
}

int NuTo::BoundaryElement1D::GetNumNodes() const
{
    return 1;
}


int NuTo::BoundaryElement1D::GetBoundaryNodeIndex(int rBoundaryNodeIndex) const
{
    Eigen::VectorXi surfaceNodeIndices = mBaseElement->GetInterpolationType()->GetSurfaceNodeIndices(mSurfaceId);
    assert(surfaceNodeIndices.rows() == 1);

    return surfaceNodeIndices(0);
}
