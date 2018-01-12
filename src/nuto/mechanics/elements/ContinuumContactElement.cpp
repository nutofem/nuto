#include "nuto/base/ErrorEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "nuto/mechanics/elements/ContinuumContactElement.h"
#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/elements/ContinuumElementIGA.h"
#include "nuto/mechanics/elements/ElementDataBase.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/ElementOutputBase.h"
#include "nuto/mechanics/elements/ElementOutputIpData.h"
#include "nuto/mechanics/elements/EvaluateDataContinuumBoundary.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationBase.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/elements/ElementShapeFunctions.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include <eigen3/Eigen/Dense>

template <int TDimSlave, int TDimMaster>
NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::ContinuumContactElement(
        const std::vector<std::pair<const ContinuumElement<TDimSlave>*, int>>& rElementsSlave,
        Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster>*, int>, Eigen::Dynamic, Eigen::Dynamic>&
                rElementsMaster,
        const ConstitutiveBase* rConstitutiveContactLaw, int rContactAlgorithm,
        const std::function<bool(int, int)>& IsNodeOnSurface)

    : ElementBase(rElementsSlave[0].first->GetStructure(), rElementsSlave[0].first->GetElementDataType(),
                  rElementsSlave[0].first->GetIpDataType(0), rElementsSlave[0].first->GetInterpolationType())
    , mElementsSlave(rElementsSlave)
    , mElementsMaster(rElementsMaster)
    , mNumDofs(0)
    , mNumSlaveDofs(0)
    , mNumMasterDofs(0)
    , mNumSlaveNodes(0)
    , mConstitutiveContactLaw(rConstitutiveContactLaw)
    , mContactType(rContactAlgorithm)
    , mIsNodeOnSurface(IsNodeOnSurface)
{
    //    if(TDimMaster != ((rSurfaceId == -1) ? TDimSlave : TDimSlave - 1))
    //        throw MechanicsException(__PRETTY_FUNCTION__, "The dimension of master side interpolation is not
    //        correct.");

    std::cout << "Number slave elements: " << rElementsSlave.size() << std::endl;

    mDofMappingComputed = false;
    mSlaveNodesMappingComputed = false;
    mMasterNodesMappingComputed = false;

    mGapMatrix.setZero(0, 0);
    mGapMatrixMeshTying.setZero(0, 0);
    mGapMatrixPenalty.setZero(0, 0);
    mSlaveShapeFunctionsWeight.setZero(0);

    mMortarGlobalGapVector.setZero(0);
    mGlobalNodalPressure.setZero(0);

    // initialize the knot values in x direction
    Eigen::VectorXd knotsX(mElementsMaster.cols() + 1);
    int i = 0;
    for (; i < mElementsMaster.cols(); i++)
        knotsX(i) = mElementsMaster(0, i).first->GetKnots()(0, 0);

    knotsX(i) = mElementsMaster(0, i - 1).first->GetKnots()(0, 1);

    mKnots.push_back(knotsX);

    // initialize the knot values in y direction
    if (TDimSlave == 3)
    {
        Eigen::VectorXd knotsY(mElementsMaster.rows() + 1);
        int i = 0;
        for (; i < mElementsMaster.rows(); i++)
            knotsY(i) = mElementsMaster(i, 0).first->GetKnots()(1, 0);

        knotsY(i) = mElementsMaster(i - 1, 0).first->GetKnots()(1, 1);
        mKnots.push_back(knotsY);
    }

    std::cout << "KnotsX: " << mKnots[0].transpose() << std::endl;
    std::cout << "KnotsY: " << mKnots[1].transpose() << std::endl;

    for (unsigned int dim = 0; dim < mKnots.size(); dim++)
    {
        for (int i = 1; i < mKnots[dim].rows(); i++)
            if (mKnots[dim](i - 1) > mKnots[dim](i))
                throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + "Knots of inconsistent ordering.");
    }
}

template <int TDimSlave, int TDimMaster>
NuTo::Element::eElementType NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GetEnumType() const
{
    return Element::eElementType::CONTINUUMCONTACTELEMENT;
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateGlobalRowDofs(
        BlockFullVector<int>& rGlobalRowDofs) const
{
    FullVector<int, Eigen::Dynamic>& dofWiseGlobalRowDofs = rGlobalRowDofs[Node::eDof::DISPLACEMENTS];

    dofWiseGlobalRowDofs.setZero(mMappingGlobal2LocalDof.size());

    // add master dofs
    for (auto& itDofs : mMappingGlobal2LocalDof)
    {
        // itMasterDofs.second = local numbering
        // itMasterDofs.first = global numbering
        // both values are unique !!, see constructor
        dofWiseGlobalRowDofs[itDofs.second] = itDofs.first;
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::ExtractAllNecessaryDofValues(
        EvaluateDataContinuumBoundary<TDimSlave>& data,
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId)
{
    data.mNodalValues[Node::eDof::DISPLACEMENTS] =
            rElementAndSurfaceId.first->ExtractNodeValues(0, Node::eDof::DISPLACEMENTS);
    data.mNodalValues[Node::eDof::COORDINATES] =
            rElementAndSurfaceId.first->ExtractNodeValues(0, Node::eDof::COORDINATES);
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateGlobalColumnDofs(
        BlockFullVector<int>& rGlobalDofMapping) const
{
    if (this->GetNumNonlocalElements() == 0)
        CalculateGlobalRowDofs(rGlobalDofMapping);
    else
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for nonlocal elements.");
}

template <int TDimSlave, int TDimMaster>
NuTo::eError NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::Evaluate(
        const ConstitutiveInputMap& rInput,
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput)
{
    (void)rInput;

    if (mDofMappingComputed == false)
        FillMappingGlobalLocalDofs();
    if (mSlaveNodesMappingComputed == false)
        FillMappingGlobalLocalSlaveNodes();

    GetConstitutiveOutputMap(rElementOutput);

    mGapMatrix.setZero(mNumDofs, mNumSlaveNodes);
    mGapMatrixAssembled = false;

    mMortarGlobalGapVector.setZero(mNumSlaveNodes);
    mMortarGlobalGapVectorAssembled = false;

    mSlaveShapeFunctionsWeight.setZero(mNumSlaveNodes);
    mSlaveShapeFunctionsWeightAssembled = false;

    mGlobalNodalPressure.setZero(mNumSlaveNodes);
    mGlobalNodalPressureAssembled = false;

    mGapMatrixPenalty.setZero(mNumDofs, mNumSlaveNodes);
    mGapMatrixPenaltyAssembled = false;

    for (const auto& it : mElementsSlave)
    {
        EvaluateDataContinuumBoundary<TDimSlave> data; //!!REFACTOR
        ExtractAllNecessaryDofValues(data, it); //!!REFACTOR
        const InterpolationBase& interpolationType = it.first->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
        int numNodesSlave = interpolationType.GetNumNodes();
        int numDofsSlave = interpolationType.GetNumDofs();

        data.mMortarGapMatrix.setZero(numDofsSlave + mNumMasterDofs, numNodesSlave);
        data.mMortarGapMatrixPenalty.setZero(numDofsSlave + mNumMasterDofs, numNodesSlave);
        data.mMortarNodalPressure.setZero(numNodesSlave);
        data.mMortarGapVector.setZero(numNodesSlave);

        data.mJacobianByWeight.setZero(this->GetNumIntegrationPoints());
        data.mShapeFunctionsIntegral.setZero(numNodesSlave);

        int numIP = this->GetNumIntegrationPoints(); // REFACTOR - the number of integration points must be at least
        // sufficient for the lower integration

        for (int theIP = 0; theIP < numIP; theIP++)
        {
            Eigen::VectorXd coordinatesIPSlave;
            Eigen::VectorXd ipCoordsNaturalSlave;
            GetGlobalIntegrationPointCoordinatesAndParameters(theIP, coordinatesIPSlave, ipCoordsNaturalSlave, it);

            auto jacobianSurface = it.first->CalculateJacobianSurface(
                    ipCoordsNaturalSlave, it.first->ExtractNodeValues(0, Node::eDof::COORDINATES), it.second);
            data.mJacobianByWeight(theIP) = jacobianSurface.norm() * this->GetIntegrationPointWeight(theIP);

            Eigen::VectorXd shapeFunsSlave = interpolationType.CalculateShapeFunctions(ipCoordsNaturalSlave);

            data.mShapeFunctionsIntegral += shapeFunsSlave * data.mJacobianByWeight(theIP);
        }

        // nodes of one element
        //        for (int k = 0; k < data.mShapeFunctionsIntegral.rows(); k++)
        //        {
        //            int surfaceid = it.second;
        //            if (!mIsNodeOnSurface(k, surfaceid))
        //                data.mShapeFunctionsIntegral(k) = 0.;
        //        }

        //        mWeightsVectorTemp = data.mShapeFunctionsIntegral;

        bool internalGradient = false;
        bool hessian_0_time_derivative = false;

        for (auto it : rElementOutput)
        {
            switch (it.first)
            {
            case Element::eOutput::INTERNAL_GRADIENT:
            {
                for (auto dofRow : this->mInterpolationType->GetActiveDofs())
                {
                    switch (dofRow)
                    {
                    case Node::eDof::DISPLACEMENTS:
                    {
                        internalGradient = true;
                        break;
                    }
                    default:
                        throw MechanicsException(__PRETTY_FUNCTION__, "Element output CONTACT_FORCE for " +
                                                                              Node::DofToString(dofRow) +
                                                                              " not implemented.");
                    }
                }
                break;
            }
            case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
            {
                for (auto dofRow : this->mInterpolationType->GetActiveDofs())
                {
                    for (auto dofCol : this->mInterpolationType->GetActiveDofs())
                    {
                        switch (Node::CombineDofs(dofRow, dofCol))
                        {
                        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
                        {
                            hessian_0_time_derivative = true;
                            break;
                        }
                        default:
                            throw MechanicsException(__PRETTY_FUNCTION__,
                                                     "Element output CONTACT_FORCE_DERIVATIVE for " +
                                                             Node::DofToString(dofRow) + " not implemented.");
                        }
                    }
                }
            }
            break;
            case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap() ...
            case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap() ...
            case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE: // already set to zero in
            // GetConstitutiveOutputMap() ...
            case Element::eOutput::GLOBAL_ROW_DOF:
            case Element::eOutput::GLOBAL_COLUMN_DOF:
            case Element::eOutput::UPDATE_STATIC_DATA:
            case Element::eOutput::UPDATE_TMP_STATIC_DATA:
            case Element::eOutput::IP_DATA:
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
            }
        }

        //===> calculate the gap matrix <===//
        GapMatrixMortarContact(data, it, internalGradient, hessian_0_time_derivative);
    }

    CalculateElementOutputs(rElementOutput);

    return eError::SUCCESSFUL;
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GapMatrixMortarTying(
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId, Eigen::MatrixXd& D,
        Eigen::MatrixXd& M)
{
    int numIP = this->GetNumIntegrationPoints();
    for (int rTheIP = 0; rTheIP < numIP; rTheIP++)
    {
        Eigen::VectorXd coordinatesIPSlave;
        Eigen::VectorXd ipCoordsNaturalSlave;
        GetGlobalIntegrationPointCoordinatesAndParameters(rTheIP, coordinatesIPSlave, ipCoordsNaturalSlave,
                                                          rElementAndSurfaceId);

        auto jacobianSurface = rElementAndSurfaceId.first->CalculateJacobianSurface(
                ipCoordsNaturalSlave, rElementAndSurfaceId.first->ExtractNodeValues(0, Node::eDof::COORDINATES),
                rElementAndSurfaceId.second);

        double jacobianByWeight = jacobianSurface.norm() * this->GetIntegrationPointWeight(rTheIP);

        // ***  Projection of the rTheIP on the master element => \xi^s_{IP}, \xi^m_*, n^m_*  *** //

        Eigen::VectorXd r;
        Eigen::VectorXd parameterMinMaster;
        const ContinuumElementIGA<TDimMaster>* masterElement =
                Projection(coordinatesIPSlave, rElementAndSurfaceId, r, parameterMinMaster);

        // *** indices for assembly - element wise *** //

        Eigen::VectorXi indicesNodesSlave, indicesNodesMaster;
        ComputeIndicesForElementAssemblyMeshTying(rElementAndSurfaceId, masterElement, indicesNodesSlave,
                                                  indicesNodesMaster);

        // *** Build the gap matrices for slave and master and assemble into D and M *** //

        const InterpolationBase& interpolationTypeDispSlave =
                rElementAndSurfaceId.first->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
        Eigen::VectorXd shapeFunsSlave = interpolationTypeDispSlave.CalculateShapeFunctions(ipCoordsNaturalSlave);

        const InterpolationBase& interpolationTypeDispMaster =
                masterElement->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
        Eigen::VectorXd shapeFunsMaster =
                interpolationTypeDispMaster.CalculateShapeFunctions(parameterMinMaster); // master is always iga

        int numSlaveFuns = shapeFunsSlave.rows();
        int numMasterFuns = shapeFunsMaster.rows();

        Eigen::MatrixXd NContactSlave;
        NContactSlave.resize(numSlaveFuns, numSlaveFuns);

        Eigen::MatrixXd NContactMaster;
        NContactMaster.resize(numSlaveFuns, numMasterFuns);

        NContactSlave = jacobianByWeight * shapeFunsSlave * shapeFunsSlave.transpose();
        NContactMaster = -jacobianByWeight * shapeFunsSlave * shapeFunsMaster.transpose();

        // *** assembly *** //
        for (int i = 0; i < indicesNodesSlave.rows(); i++)
            for (int j = 0; j < indicesNodesSlave.rows(); j++)
                D(indicesNodesSlave(i), indicesNodesSlave(j)) += NContactSlave(i, j);

        for (int i = 0; i < indicesNodesSlave.rows(); i++)
            for (int j = 0; j < indicesNodesMaster.rows(); j++)
                M(indicesNodesSlave(i), indicesNodesMaster(j)) += NContactMaster(i, j);

        for (int i = 0; i < indicesNodesSlave.rows(); i++)
            mSlaveShapeFunctionsWeight(indicesNodesSlave(i)) += shapeFunsSlave(i) * jacobianByWeight;
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GapMatrixMortarContact(
        EvaluateDataContinuumBoundary<TDimSlave>& rData,
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId, bool rAssembleForce,
        bool rAssembleDerivative)
{
    int numIP = this->GetNumIntegrationPoints();
    for (int rTheIP = 0; rTheIP < numIP; rTheIP++)
    {
        Eigen::VectorXd coordinatesIPSlave;
        Eigen::VectorXd ipCoordsNaturalSlave;
        GetGlobalIntegrationPointCoordinatesAndParameters(rTheIP, coordinatesIPSlave, ipCoordsNaturalSlave,
                                                          rElementAndSurfaceId);

        auto jacobianSurface = rElementAndSurfaceId.first->CalculateJacobianSurface(
                ipCoordsNaturalSlave, rElementAndSurfaceId.first->ExtractNodeValues(0, Node::eDof::COORDINATES),
                rElementAndSurfaceId.second);

        double jacobianByWeight = jacobianSurface.norm() * this->GetIntegrationPointWeight(rTheIP);

        // ***  Projection of the rTheIP on the master element => \xi^s_{IP}, \xi^m_*, n^m_*  *** //

        Eigen::VectorXd r;
        Eigen::VectorXd parameterMinMaster;
        const ContinuumElementIGA<TDimMaster>* masterElement =
                Projection(coordinatesIPSlave, rElementAndSurfaceId, r, parameterMinMaster);

        // *** Build the gap matrix *** //

        // normal vector
        Eigen::VectorXd normal = masterElement->InterpolateDofGlobalSurfaceNormal(parameterMinMaster);

        // Assemble the element gap matrix => \int_{\Gamma_e} F(ShapeFunctionsSlave(\xi^s),
        // ShapeFunctionsMaster(\xi^*),
        // n^*) d\Gamma

        const InterpolationBase& interpolationTypeDispSlave =
                rElementAndSurfaceId.first->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
        Eigen::VectorXd shapeFunsSlave = interpolationTypeDispSlave.CalculateShapeFunctions(ipCoordsNaturalSlave);

        const InterpolationBase& interpolationTypeDispMaster =
                masterElement->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
        Eigen::VectorXd shapeFunsMaster =
                interpolationTypeDispMaster.CalculateShapeFunctions(parameterMinMaster); // master is always iga


        int numSlaveFunsNormal = shapeFunsSlave.rows() * normal.rows();
        int numMasterFunsNormal = shapeFunsMaster.rows() * normal.rows();

        Eigen::VectorXd shapeFunsSlaveNormal(numSlaveFunsNormal);
        int count = 0;
        for (int i = 0; i < shapeFunsSlave.rows(); i++)
        {
            shapeFunsSlaveNormal.block(count, 0, normal.rows(), 1) = shapeFunsSlave(i) * normal;
            count += normal.rows();
        }

        Eigen::VectorXd shapeFunsMasterNormal(numMasterFunsNormal);
        count = 0;
        for (int i = 0; i < shapeFunsMaster.rows(); i++)
        {
            shapeFunsMasterNormal.block(count, 0, normal.rows(), 1) = (-1) * shapeFunsMaster(i) * normal;
            count += normal.rows();
        }

        double gap = r.dot(normal); // gap at ip
        //        std::cout << rTheIP << ": " << gap << std::endl << std::flush;

        // ----------- geometrical linearization terms ---- start ------------//

        if (0)
        {
            Eigen::VectorXd shapeFunsMasterFirstDer =
                    interpolationTypeDispMaster.CalculateDerivativeShapeFunctionsNatural(
                            parameterMinMaster); // master is always iga

            Eigen::VectorXd tau_1 = masterElement->InterpolateDofGlobalSurfaceDerivative(0, parameterMinMaster, 1, 0);
            double m_11 = tau_1.dot(tau_1);
            Eigen::VectorXd secondder_1 =
                    masterElement->InterpolateDofGlobalSurfaceDerivative(0, parameterMinMaster, 2, 0);
            double k_11 = secondder_1.dot(normal);

            Eigen::VectorXd T(numSlaveFunsNormal + numMasterFunsNormal);
            Eigen::VectorXd N_1(numSlaveFunsNormal + numMasterFunsNormal);

            T.setZero();
            N_1.setZero();

            count = 0;

            // slave part
            for (int i = 0; i < shapeFunsSlave.rows(); i++)
            {
                T.block(count, 0, tau_1.rows(), 1) = shapeFunsSlave(i) * tau_1;
                count += tau_1.rows();
            }

            // master part
            for (int i = 0; i < shapeFunsSlave.rows(); i++)
            {
                T.block(count, 0, tau_1.rows(), 1) = shapeFunsMaster(i) * tau_1;
                N_1.block(count, 0, tau_1.rows(), 1) = shapeFunsMasterFirstDer(i) * normal;
                count += tau_1.rows();
            }

            Eigen::VectorXd D_1 = (1 / (m_11 - k_11 * gap)) * (T - gap * N_1);
            Eigen::VectorXd N_1_square = N_1 - k_11 * D_1;

            Eigen::MatrixXd L = (-gap / m_11) * N_1_square * N_1_square.transpose() - D_1 * N_1.transpose() -
                                N_1 * D_1.transpose() + k_11 * D_1 * D_1.transpose();
        }

        // ----------- geometrical linearization terms --- end -------------//

        // matrix containing the derivatives of the slave and master side multiplied by the normal (R^S * normal
        // \\  R^M
        // * normal)
        Eigen::MatrixXd NContact;
        NContact.resize(numSlaveFunsNormal + numMasterFunsNormal, shapeFunsSlave.rows());

        NContact.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows()) =
                shapeFunsSlaveNormal * shapeFunsSlave.transpose();
        NContact.block(numSlaveFunsNormal, 0, numMasterFunsNormal, shapeFunsSlave.rows()) =
                shapeFunsMasterNormal * shapeFunsSlave.transpose();
        NContact *= jacobianByWeight;

        double derivativeContactForce = mConstitutiveContactLaw->GetContactForceDerivative(gap);

        // GAP-MATRIX
        rData.mMortarGapMatrix.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows()) +=
                NContact.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows());

        if (mContactType == 1 && rAssembleDerivative)
            rData.mMortarGapMatrixPenalty.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows()) +=
                    NContact.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows()) * derivativeContactForce;

        const int numNodes = interpolationTypeDispMaster.GetNumNodes();
        unsigned int numDofsPerType =
                masterElement->GetNode(interpolationTypeDispMaster.GetNodeIndex(0))->GetNum(Node::eDof::DISPLACEMENTS);
        count = 0;
        for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
        {
            const NodeBase* nodePtr = masterElement->GetNode(interpolationTypeDispMaster.GetNodeIndex(iNodeDof));

            for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
            {
                int dofID = nodePtr->GetDof(Node::eDof::DISPLACEMENTS, iDof);
                int index = mMappingGlobal2LocalDof[dofID];
                rData.mMortarGapMatrix.row(numSlaveFunsNormal + index - mNumSlaveDofs) +=
                        NContact.row(count + numSlaveFunsNormal);

                if (mContactType == 1 && rAssembleDerivative)
                    rData.mMortarGapMatrixPenalty.row(numSlaveFunsNormal + index - mNumSlaveDofs) +=
                            NContact.row(count + numSlaveFunsNormal) * derivativeContactForce;

                count++;
            }
        }

        // *** penalty vector or gap vector for the contact force *** //
        if (mContactType == 0)
        {
            // GAP-VECTOR
            rData.mMortarGapVector += shapeFunsSlave * gap * jacobianByWeight;
        }
        else if (mContactType == 1 && rAssembleForce)
        {
            // NODAL-PRESSURE
            rData.mMortarNodalPressure +=
                    shapeFunsSlave * (mConstitutiveContactLaw->GetContactForce(gap)) * jacobianByWeight;
        }
    }

    //    mGapMatrixTemp =
    //            rData.mMortarGapMatrix.block(0, 0, rData.mMortarGapVector.rows() * 3, rData.mMortarGapVector.rows());

    //    mGapVectorTemp = rData.mMortarGapVector;

    if (mContactType == 0)
    {
        AssembleMortarTypeContact(rData, rElementAndSurfaceId);
    }
    else if (mContactType == 1)
    {
        AssembleNonMortarTypeContact(rData, rElementAndSurfaceId, rAssembleForce, rAssembleDerivative);
    }
}

template <int TDimSlave, int TDimMaster>
const NuTo::ContinuumElementIGA<TDimMaster>* NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::Projection(
        const Eigen::VectorXd& coordinatesIPSlave,
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId,
        Eigen::VectorXd& rProjectionVector, Eigen::VectorXd& rParameterMinMaster)
{
    // **** Get the starting point for iteration **** //

    int numStartingPointsOneDir = 10;
    Eigen::VectorXd coords(numStartingPointsOneDir);
    coords(0) = -0.99;
    double length = 1.98;
    for (int i = 1; i < numStartingPointsOneDir; i++)
        coords(i) = coords(i - 1) + length / (numStartingPointsOneDir - 1);

    int numStartingPointsElement = std::pow(numStartingPointsOneDir, TDimMaster);
    std::vector<Eigen::Matrix<double, TDimMaster, 1>> referenceCoordinates(numStartingPointsElement);

    if (mKnots.size() == 1)
    {
        int count = 0;
        for (int x = 0; x < numStartingPointsOneDir; x++)
        {
            referenceCoordinates[count](0, 0) = coords(x);
            count++;
        }
    }

    if (mKnots.size() == 2)
    {
        int count = 0;

        for (int x = 0; x < numStartingPointsOneDir; x++)
        {
            for (int y = 0; y < numStartingPointsOneDir; y++)
            {
                referenceCoordinates[count](0, 0) = coords(x);
                referenceCoordinates[count](1, 0) = coords(y);
                count++;
            }
        }
    }

    double minDistance = std::numeric_limits<double>::infinity();
    Eigen::Vector2i indexMasterElement(0, 0);
    Eigen::Matrix<double, TDimMaster, 1> localParameter;
    for (int i = 0; i < mElementsMaster.rows(); i++) // y
    {
        for (int j = 0; j < mElementsMaster.cols(); j++) // x
        {
            const auto* elementPtr = mElementsMaster(i, j).first;
            int surfaceId = mElementsMaster(i, j).second;

            // ===> Get the position on the master curve/surface
            // ===> Compare and set to minimum if though
            const InterpolationBase& interpolationTypeDispMaster =
                    elementPtr->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
            for (auto it : referenceCoordinates)
            {
                Eigen::VectorXd parameter = interpolationTypeDispMaster.CalculateNaturalSurfaceCoordinates(
                        it, surfaceId, elementPtr->GetKnots());
                Eigen::VectorXd coordinatesMaster =
                        elementPtr->InterpolateDofGlobalSurfaceDerivative(0, parameter, 0, 0);

                double distance = (coordinatesMaster - coordinatesIPSlave).norm();
                if (minDistance > distance)
                {
                    localParameter = it;
                    minDistance = distance;
                    rParameterMinMaster = parameter;
                    indexMasterElement(0) = i; // y
                    indexMasterElement(1) = j; // x
                }
            }
        }
    }

    const NuTo::ContinuumElementIGA<TDimMaster>* masterElement =
            mElementsMaster(indexMasterElement(0), indexMasterElement(1)).first;

    int masterSurfaceId = mElementsMaster(indexMasterElement(0), indexMasterElement(1)).second;

    // **** Newton iteration ****//
    double tol = 1.e-8;
    double error = 1.;
    int maxNumIter = 100;
    int numIter = 0;
    Eigen::VectorXd coordinatesMaster;

    int numPrimes = (rElementAndSurfaceId.second == -1) ? TDimSlave : TDimSlave - 1;

    while (error > tol && numIter < maxNumIter)
    {
        // ==> function (dprime)
        coordinatesMaster = masterElement->InterpolateDofGlobalSurfaceDerivative(0, rParameterMinMaster, 0, 0);
        rProjectionVector = coordinatesIPSlave - coordinatesMaster;
        Eigen::VectorXd dprime(numPrimes, 1);
        dprime.setZero(numPrimes);
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> prime =
                masterElement->InterpolateDofGlobalSurfaceDerivativeTotal(0, rParameterMinMaster, 1, masterSurfaceId);
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> primeprime =
                masterElement->InterpolateDofGlobalSurfaceDerivativeTotal(0, rParameterMinMaster, 2, masterSurfaceId);

        for (int j = 0; j < prime.cols(); j++)
            dprime(j) += rProjectionVector.dot(prime(0, j));

        // ==> derivative
        Eigen::MatrixXd dprimeprime(numPrimes, numPrimes);
        dprimeprime.setZero(numPrimes, numPrimes);

        for (int i = 0; i < prime.cols(); i++)
        {
            for (int j = 0; j < prime.cols(); j++)
            {
                dprimeprime(i, j) = rProjectionVector.dot(primeprime(i, j)) - prime(0, i).dot(prime(0, j));
            }
        }

        // iteration step
        Eigen::VectorXd increment = dprimeprime.colPivHouseholderQr().solve(-dprime);

        // New parameter value (check in which element)
        // for the FindSpan function degree = 0, since no multiple knots at the beginning and end
        if (mKnots.size() == 1)
        {
            rParameterMinMaster(0) += increment(0);
            if (rParameterMinMaster(0) < 0.)
                rParameterMinMaster(0) = 0.;

            if (rParameterMinMaster(0) > 1.)
                rParameterMinMaster(0) = 1.;
            indexMasterElement(1) = ShapeFunctionsIGA::FindSpan(rParameterMinMaster(0), 0, mKnots[0]);
        }
        else if (mKnots.size() == 2)
        {
            //            if (rParameterMinMaster(0) + increment(0) < 0. || rParameterMinMaster(1) + increment(1) < 0.
            //            ||
            //                rParameterMinMaster(0) + increment(0) > 1. || rParameterMinMaster(1) + increment(1) > 1.)
            //            {
            //                std::cout << "dprime.norm: " << dprime.norm() << ", increment.norm: " << increment.norm()
            //                << std::endl;
            //                break;
            //            }

            rParameterMinMaster += increment;

            indexMasterElement(1) = ShapeFunctionsIGA::FindSpan(rParameterMinMaster(0), 0, mKnots[0]); // x
            indexMasterElement(0) = ShapeFunctionsIGA::FindSpan(rParameterMinMaster(1), 0, mKnots[1]); // y
        }

        masterElement = mElementsMaster(indexMasterElement(0), indexMasterElement(1)).first;

        error = dprime.norm();
        numIter++;
    }

    //    if (numIter >= maxNumIter)
    //        std::cout << "!!!!!!ContinuumContactElement: Maximum number of Newton iterations exceeded! (error = " <<
    //        error
    //                  << ")" << std::endl;

    Eigen::VectorXd normal = masterElement->InterpolateDofGlobalSurfaceNormal(rParameterMinMaster);
    Eigen::VectorXd point = masterElement->InterpolateDofGlobalSurfaceDerivative(0, rParameterMinMaster, 0, 0);
    Eigen::VectorXd slavePointIteration = point + (rProjectionVector.dot(normal)) * normal;

    if ((slavePointIteration - coordinatesIPSlave).lpNorm<Eigen::Infinity>() > 1.e-8)
    {
        std::cout << "Not converged, deviation L2 norm (" << coordinatesIPSlave.transpose()
                  << "): " << (slavePointIteration - coordinatesIPSlave).norm() << std::endl;
    }

    return masterElement;
}

// Assembling: A_A, g_{NA}, G
template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::AssembleMortarTypeContact(
        const EvaluateDataContinuumBoundary<TDimSlave>& rData,
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId, bool rAssembleGapVector)
{
    // assemble into the contact force
    // ==> slave
    auto dof = Node::eDof::DISPLACEMENTS;

    const InterpolationType* interpolationType = rElementAndSurfaceId.first->GetInterpolationType();
    const InterpolationBase& interpolationDof = interpolationType->Get(dof);
    const int numNodes = interpolationDof.GetNumNodes();

    unsigned int numDofsPerType = rElementAndSurfaceId.first->GetNode(interpolationDof.GetNodeIndex(0))->GetNum(dof);

    Eigen::VectorXi indicesDofs;
    indicesDofs.setZero(numNodes * numDofsPerType + mNumMasterDofs);

    Eigen::VectorXi indicesNodes;
    indicesNodes.setZero(numNodes);

    // *** penalty vector or gap vector for the contact force ***
    int countDofs = 0;
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const NodeBase* nodePtr = rElementAndSurfaceId.first->GetNode(interpolationDof.GetNodeIndex(iNode));

        int globalNodeId = this->mStructure->NodeGetId(nodePtr);
        int indexNode = mMappingGlobal2LocalSlaveNodes[globalNodeId];
        indicesNodes(iNode) = indexNode;

        mSlaveShapeFunctionsWeight(indexNode) += rData.mShapeFunctionsIntegral(iNode);

        if (rAssembleGapVector)
            mMortarGlobalGapVector(indexNode) += rData.mMortarGapVector(iNode);

        for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
        {
            int indexDofs = mMappingGlobal2LocalDof[nodePtr->GetDof(Node::eDof::DISPLACEMENTS, iDof)];
            indicesDofs(countDofs) = indexDofs;
            countDofs++;
        }
    }

    int countMaster = mNumSlaveDofs;
    // add master dofs (all master dofs in every contact force vector/matrix)
    for (int i = 0; i < mNumMasterDofs; i++)
    {
        indicesDofs(countDofs) = countMaster;
        countDofs++;
        countMaster++;
    }

    for (int i = 0; i < indicesDofs.rows(); i++)
    {
        for (int j = 0; j < indicesNodes.rows(); j++)
        {
            mGapMatrix(indicesDofs(i), indicesNodes(j)) += rData.mMortarGapMatrix(i, j);
        }
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::AssembleNonMortarTypeContact(
        const EvaluateDataContinuumBoundary<TDimSlave>& rData,
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId, bool rAssembleForce,
        bool rAssembleDerivative)
{
    // assemble into the contact force
    // ==> slave
    auto dof = Node::eDof::DISPLACEMENTS;

    const InterpolationType* interpolationType = rElementAndSurfaceId.first->GetInterpolationType();
    const InterpolationBase& interpolationDof = interpolationType->Get(dof);
    const int numNodes = interpolationDof.GetNumNodes();

    unsigned int numDofsPerType = rElementAndSurfaceId.first->GetNode(interpolationDof.GetNodeIndex(0))->GetNum(dof);

    Eigen::VectorXi indicesDofs;
    indicesDofs.setZero(numNodes * numDofsPerType + mNumMasterDofs);

    Eigen::VectorXi indicesNodes;
    indicesNodes.setZero(numNodes);

    // *** penalty vector or gap vector for the contact force ***
    int countNode = 0;
    int countDofs = 0;
    for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
    {
        const NodeBase* nodePtr = rElementAndSurfaceId.first->GetNode(interpolationDof.GetNodeIndex(iNodeDof));

        int globalNodeId = this->mStructure->NodeGetId(nodePtr);
        int indexNode = mMappingGlobal2LocalSlaveNodes[globalNodeId];
        indicesNodes(countNode) = indexNode;

        mSlaveShapeFunctionsWeight(indexNode) += rData.mShapeFunctionsIntegral(countNode);

        if (rAssembleForce)
            mGlobalNodalPressure(indexNode) += rData.mMortarNodalPressure(countNode);

        countNode++;

        for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
        {
            int indexDofs = mMappingGlobal2LocalDof[nodePtr->GetDof(Node::eDof::DISPLACEMENTS, iDof)];
            indicesDofs(countDofs) = indexDofs;
            countDofs++;
        }
    }

    int countMaster = mNumSlaveDofs;
    // add master dofs (all master dofs in every contact force vector/matrix)
    for (int i = 0; i < mNumMasterDofs; i++)
    {
        indicesDofs(countDofs) = countMaster;
        countDofs++;
        countMaster++;
    }

    for (int i = 0; i < indicesDofs.rows(); i++)
    {
        for (int j = 0; j < indicesNodes.rows(); j++)
        {
            mGapMatrix(indicesDofs(i), indicesNodes(j)) += rData.mMortarGapMatrix(i, j);
            if (rAssembleDerivative)
                mGapMatrixPenalty(indicesDofs(i), indicesNodes(j)) += rData.mMortarGapMatrixPenalty(i, j);
        }
    }
}

template <int TDimSlave, int TDimMaster>
NuTo::ConstitutiveOutputMap NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GetConstitutiveOutputMap(
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) const
{
    ConstitutiveOutputMap constitutiveOutput;

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::eOutput::INTERNAL_GRADIENT:
        case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
            // no resize needed, since the mortar matrices are set at the end of the evaluate routine...
            break;
        case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE: // not needed, so set to 0
            FillConstitutiveOutputMapHessian1(constitutiveOutput, it.second->GetBlockFullMatrixDouble());
            break;
        case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE: // not needed, so set to 0
            FillConstitutiveOutputMapHessian2Lumped(constitutiveOutput, it.second->GetBlockFullVectorDouble());
            break;
        case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE: // not needed, so set to 0
            FillConstitutiveOutputMapHessian2(constitutiveOutput, it.second->GetBlockFullMatrixDouble());
            break;
        case Element::eOutput::UPDATE_STATIC_DATA:
            constitutiveOutput[Constitutive::eOutput::UPDATE_STATIC_DATA] = 0;
            break;
        case Element::eOutput::UPDATE_TMP_STATIC_DATA:
            constitutiveOutput[Constitutive::eOutput::UPDATE_TMP_STATIC_DATA] = 0;
            break;
        case Element::eOutput::IP_DATA:
            //            this->FillConstitutiveOutputMapIpData(constitutiveOutput, it.second->GetIpData());
            break;
        case Element::eOutput::GLOBAL_ROW_DOF:
            CalculateGlobalRowDofs(it.second->GetBlockFullVectorInt());
            break;
        case Element::eOutput::GLOBAL_COLUMN_DOF:
            CalculateGlobalColumnDofs(it.second->GetBlockFullVectorInt());
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element  output not implemented.");
        }
    }
    return constitutiveOutput;
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillConstitutiveOutputMapHessian1(
        ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian1) const
{
    (void)rConstitutiveOutput;
    auto activeDofs = this->mInterpolationType->GetActiveDofs();
    if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact is only implemented for displacements.");
    rHessian1(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS).Resize(mNumDofs, mNumDofs);
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillConstitutiveOutputMapHessian2Lumped(
        ConstitutiveOutputMap& rConstitutiveOutput, BlockFullVector<double>& rHessian2Lumped) const
{
    (void)rConstitutiveOutput;
    auto activeDofs = this->mInterpolationType->GetActiveDofs();
    if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact is only implemented for displacements.");
    rHessian2Lumped[Node::eDof::DISPLACEMENTS].Resize(mNumDofs);
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillConstitutiveOutputMapHessian2(
        ConstitutiveOutputMap& rConstitutiveOutput, BlockFullMatrix<double>& rHessian2) const
{
    (void)rConstitutiveOutput;
    auto activeDofs = this->mInterpolationType->GetActiveDofs();
    if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact is only implemented for displacements.");
    rHessian2(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS).Resize(mNumDofs, mNumDofs);
}


template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateElementOutputs(
        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput)
{
    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::eOutput::INTERNAL_GRADIENT:
            CalculateElementOutputContactForce(it.second->GetBlockFullVectorDouble());
            break;
        case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
            CalculateElementOutputContactForceDerivative(it.second->GetBlockFullMatrixDouble());
            break;
        case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
        case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
        case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE: // already set to zero in
        // GetConstitutiveOutputMap()...
        case Element::eOutput::GLOBAL_ROW_DOF:
        case Element::eOutput::GLOBAL_COLUMN_DOF:
        case Element::eOutput::UPDATE_STATIC_DATA:
        case Element::eOutput::UPDATE_TMP_STATIC_DATA:
        case Element::eOutput::IP_DATA:
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
        }
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateElementOutputContactForce(
        BlockFullVector<double>& rInternalGradient)
{
    for (auto dofRow : this->mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
        {
            for (int i = 0; i < mNumSlaveNodes; i++)
            {
                //                std::cout << mSlaveShapeFunctionsWeight(i) << " ";
                if (fabs(mSlaveShapeFunctionsWeight(i)) < 1.e-10)
                {
                    mSlaveShapeFunctionsWeight(i) = 1.;
                    //                    mWeightsVectorTemp(i) = 1.;
                }
            }

            //            std::cout << mSlaveShapeFunctionsWeight << std::endl;

            //            Eigen::MatrixXd mat = ((mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal());
            //            std::cout << mat << std::endl;

            //            std::cout << "Gap-Matri: \n" << std::endl;
            //            std::cout << mGapMatrix << std::endl;

            Eigen::MatrixXd gapMatrixScaled = mGapMatrix * ((mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal());

            //            std::cout << "Gap-Matri-Scaled: \n" << std::endl;
            //            std::cout << gapMatrixScaled << std::endl;

            if (mContactType == 0)
            {
                for (int i = 0; i < mNumSlaveNodes; i++)
                    mGlobalNodalPressure(i) = mConstitutiveContactLaw->GetContactForce(mMortarGlobalGapVector(i));

                //                mPressureVectorTemp.setZero(mNumSlaveNodes);
                //                for (int i = 0; i < mNumSlaveNodes; i++)
                //{
                //                    mPressureVectorTemp(i) =
                //                    mConstitutiveContactLaw->GetContactForce(mGapVectorTemp(i));
                //}

                //                int numForces = 0;
                //                for (int i = 0; i < mGlobalNodalPressure.rows(); i++)
                //                {
                //                    if (mGlobalNodalPressure(i) < 0.)
                //                        numForces++;
                //                }
                //                std::cout << numForces << std::endl;

                Eigen::VectorXd force = gapMatrixScaled * mGlobalNodalPressure;

                //                Eigen::VectorXd forceTemp =
                //                        mGapMatrixTemp * ((mWeightsVectorTemp.cwiseInverse()).asDiagonal()) *
                //                        mPressureVectorTemp;

                //                std::cout << (forceTemp - force.block(0, 0, forceTemp.rows(), 1)).norm() << std::endl;

                //                Eigen::VectorXd forceTemp =
                //                        (mGapMatrixTemp * ((mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal()))
                //                        *
                //                        mGlobalNodalPressure;

                rInternalGradient[dofRow] = gapMatrixScaled * mGlobalNodalPressure;

                std::cout << "|force|_1: ";
                double sum = 0.;
                for (int i = 0; i < mNumSlaveDofs; i++)
                {
                    sum += force(i);
                }
                std::cout << sum << std::endl;
            }
            else if (mContactType == 1)
            {
                rInternalGradient[dofRow] = gapMatrixScaled * mGlobalNodalPressure;
            }
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Element output CONTACT_FORCE for " +
                                                                  Node::DofToString(dofRow) + " not implemented.");
        }
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateElementOutputContactForceDerivative(
        BlockFullMatrix<double>& rGapMatrix)
{
    for (auto dofRow : this->mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : this->mInterpolationType->GetActiveDofs())
        {
            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            {
                for (int i = 0; i < mNumSlaveNodes; i++)
                    if (fabs(mSlaveShapeFunctionsWeight(i)) < 1.e-10)
                        mSlaveShapeFunctionsWeight(i) = 1.;

                Eigen::MatrixXd gapMatrixScaledWithIntegrals =
                        mGapMatrix * ((mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal());

                if (mContactType == 0)
                {
                    Eigen::VectorXd globalNodalPressureDerivative(mNumSlaveNodes);
                    for (int i = 0; i < mNumSlaveNodes; i++)
                        globalNodalPressureDerivative(i) =
                                mConstitutiveContactLaw->GetContactForceDerivative(mMortarGlobalGapVector(i));

                    Eigen::MatrixXd gapMatrixScaledWithForceDers =
                            mGapMatrix * (globalNodalPressureDerivative.asDiagonal());

                    // Eigen::MatrixXd mat =
                    // gapMatrixScaledWithIntegrals*(gapMatrixScaledWithForceDers.transpose());

                    rGapMatrix(dofRow, dofCol) =
                            gapMatrixScaledWithIntegrals * (gapMatrixScaledWithForceDers.transpose());
                    //                    rGapMatrix(dofRow, dofCol) =
                    //                    mGapMatrix*((mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal())*mGapMatrix*(globalNodalPressureDerivative.asDiagonal());
                }
                else if (mContactType == 1)
                {
                    rGapMatrix(dofRow, dofCol) = gapMatrixScaledWithIntegrals * (mGapMatrixPenalty.transpose());
                    //                    rGapMatrix(dofRow, dofCol) =
                    //                    mGapMatrix*((mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal())*(mGapMatrixPenalty.transpose());
                }
                break;
            }
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Element output CONTACT_FORCE_DERIVATIVE for " +
                                                                      Node::DofToString(dofRow) + " not implemented.");
            }
        }
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GetGlobalIntegrationPointCoordinatesAndParameters(
        int rIpNum, Eigen::VectorXd& rCoordinatesIPSlave, Eigen::VectorXd& rParamsIPSlave,
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rElementAndSurfaceId) const
{
    Eigen::VectorXd naturalSurfaceIpCoordinates;
    if (mKnots.size() == 1)
    {
        double ipCoordinate;
        this->GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, ipCoordinate);
        naturalSurfaceIpCoordinates.resize(1);
        naturalSurfaceIpCoordinates(0) = ipCoordinate;
    }
    else if (mKnots.size() == 2)
    {
        double ipCoordinates[2];
        const NuTo::IntegrationTypeBase* ptrInType = this->GetIntegrationType();
        ptrInType->GetLocalIntegrationPointCoordinates2D(rIpNum, ipCoordinates);
        naturalSurfaceIpCoordinates.resize(2);
        naturalSurfaceIpCoordinates(0) = ipCoordinates[0];
        naturalSurfaceIpCoordinates(1) = ipCoordinates[1];
    }

    // ===> Get the position \xi^s_{IP}
    const InterpolationBase& interpolationTypeCoordsSlave =
            rElementAndSurfaceId.first->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);

    if (interpolationTypeCoordsSlave.GetTypeOrder() == Interpolation::eTypeOrder::SPLINE)
        rParamsIPSlave = interpolationTypeCoordsSlave.CalculateNaturalSurfaceCoordinates(
                naturalSurfaceIpCoordinates, rElementAndSurfaceId.second,
                rElementAndSurfaceId.first->GetKnots()); // IGA
    else
        rParamsIPSlave = interpolationTypeCoordsSlave.CalculateNaturalSurfaceCoordinates(
                naturalSurfaceIpCoordinates, rElementAndSurfaceId.second); // FEM

    rCoordinatesIPSlave = rElementAndSurfaceId.first->InterpolateDofGlobalCurrentConfiguration(
            0, rParamsIPSlave, Node::eDof::COORDINATES, Node::eDof::DISPLACEMENTS);
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::ComputeIndicesForElementAssemblyMeshTying(
        const std::pair<const ContinuumElement<TDimSlave>*, int>& rSlaveElementAndSurfaceId,
        const ContinuumElement<TDimMaster>* rMasterElement, Eigen::VectorXi& indicesNodesSlave,
        Eigen::VectorXi& indicesNodesMaster)
{
    auto dof = Node::eDof::DISPLACEMENTS; // only disps for contact allowed

    // slave
    const InterpolationType* interpolationTypeSlave = rSlaveElementAndSurfaceId.first->GetInterpolationType();
    const InterpolationBase& interpolationDofSlave = interpolationTypeSlave->Get(dof);
    const int numNodesSlave = interpolationDofSlave.GetNumNodes();

    indicesNodesSlave.setZero(numNodesSlave);

    int countNodes = 0;
    for (int iNode = 0; iNode < numNodesSlave; ++iNode)
    {
        const NodeBase* nodePtrSlave =
                rSlaveElementAndSurfaceId.first->GetNode(interpolationDofSlave.GetNodeIndex(iNode));
        int globalNodeId = this->mStructure->NodeGetId(nodePtrSlave);
        int indexNode = mMappingGlobal2LocalSlaveNodes[globalNodeId];
        indicesNodesSlave(countNodes) = indexNode;
        countNodes++;
    }

    // master
    const InterpolationType* interpolationTypeMaster = rMasterElement->GetInterpolationType();
    const InterpolationBase& interpolationDofMaster = interpolationTypeMaster->Get(dof);
    const int numNodesMaster = interpolationDofMaster.GetNumNodes();

    indicesNodesMaster.setZero(numNodesMaster);

    countNodes = 0;
    for (int iNode = 0; iNode < numNodesMaster; ++iNode)
    {
        const NodeBase* nodePtrMaster = rMasterElement->GetNode(interpolationDofMaster.GetNodeIndex(iNode));
        int globalNodeId = this->mStructure->NodeGetId(nodePtrMaster);
        int indexNode = mMappingGlobal2LocalMasterNodes[globalNodeId];
        indicesNodesMaster(countNodes) = indexNode;
        countNodes++;
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::ComputeMeshTyingMatrix(
        Eigen::MatrixXd& D, Eigen::MatrixXd& M, std::unordered_map<int, int>& mappingGlobal2LocalSlaveNode,
        std::unordered_map<int, int>& mappingGlobal2LocalMasterNode)
{
    if (mSlaveNodesMappingComputed == false)
        FillMappingGlobalLocalSlaveNodes();
    if (mMasterNodesMappingComputed == false)
        FillMappingGlobalLocalMasterNodes();

    Eigen::MatrixXd D_unscaled, M_unscaled;
    D_unscaled.setZero(mNumSlaveNodes, mNumSlaveNodes);
    M_unscaled.setZero(mNumSlaveNodes, mNumMasterNodes);

    mSlaveShapeFunctionsWeight.setZero(mNumSlaveNodes);

    for (const auto& it : mElementsSlave)
    {
        EvaluateDataContinuumBoundary<TDimSlave> data; //!!REFACTOR
        ExtractAllNecessaryDofValues(data, it); //!!REFACTOR

        // *** calculate gap mortrar mesh tying matrices *** //
        GapMatrixMortarTying(it, D_unscaled, M_unscaled);
    }

    std::vector<int> nonZeroEntries;
    // *** scaling with ansatz functions *** //
    for (int i = 0; i < mNumSlaveNodes; i++)
    {
        if (fabs(mSlaveShapeFunctionsWeight(i)) < 1.e-10)
            mSlaveShapeFunctionsWeight(i) = 1.;
        else
            nonZeroEntries.push_back(i);
    }

    for (int it : nonZeroEntries)
    {
        for (auto& itMap : mMappingGlobal2LocalSlaveNodes)
        {
            if (itMap.second == it)
                mappingGlobal2LocalSlaveNode.insert(std::make_pair(itMap.first, itMap.second));
        }
    }

    mappingGlobal2LocalMasterNode = mMappingGlobal2LocalMasterNodes;

    D = D_unscaled;
    M = M_unscaled;

    // scaling - not needed becaise Variational contribution = 0, contrary to the contact
    //    D = (mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal() * D_unscaled;
    //    M = (mSlaveShapeFunctionsWeight.cwiseInverse()).asDiagonal() * M_unscaled;
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillMappingGlobalLocalDofs()
{
    mNumSlaveDofs = 0;
    for (auto& it : mElementsSlave)
    {
        auto activeDofs = it.first->GetInterpolationType()->GetActiveDofs();
        if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
            throw MechanicsException(__PRETTY_FUNCTION__, "Contact element is only implemented for displacements.");

        auto dof = Node::eDof::DISPLACEMENTS;

        const InterpolationType* interpolationType = it.first->GetInterpolationType();
        const InterpolationBase& interpolationDof = interpolationType->Get(dof);
        const int numNodes = interpolationDof.GetNumNodes();

        unsigned int numDofsPerType = it.first->GetNode(interpolationDof.GetNodeIndex(0))->GetNum(dof);

        for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
        {
            const NodeBase* nodePtr = it.first->GetNode(interpolationDof.GetNodeIndex(iNodeDof));

            for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
            {
                int dofID = nodePtr->GetDof(dof, iDof);
                std::unordered_map<int, int>::const_iterator got = mMappingGlobal2LocalDof.find(dofID);
                if (got == mMappingGlobal2LocalDof.end()) // not found => insert a new one
                {
                    mMappingGlobal2LocalDof.insert(std::make_pair(dofID, mNumSlaveDofs));
                    //                    std::cout << "(dofID, mNumSlaveDofs): " << dofID << " " << mNumSlaveDofs <<
                    //                    std::endl;
                    mNumSlaveDofs++;
                }
            }
        }
    }

    mNumMasterDofs = 0;
    int localDofIndex = mNumSlaveDofs;

    for (int i = 0; i < mElementsMaster.rows(); i++)
    {
        for (int j = 0; j < mElementsMaster.cols(); j++)
        {
            const ContinuumElementIGA<TDimMaster>* masterElement = mElementsMaster(i, j).first;
            auto activeDofs = masterElement->GetInterpolationType()->GetActiveDofs();
            if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
                throw MechanicsException(__PRETTY_FUNCTION__, "Contact element is only implemented for displacements.");

            auto dof = Node::eDof::DISPLACEMENTS;

            const InterpolationType* interpolationType = masterElement->GetInterpolationType();
            const InterpolationBase& interpolationDof = interpolationType->Get(dof);
            const int numNodes = interpolationDof.GetNumNodes();

            unsigned int numDofsPerType = masterElement->GetNode(interpolationDof.GetNodeIndex(0))->GetNum(dof);

            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                const NodeBase* nodePtr = masterElement->GetNode(interpolationDof.GetNodeIndex(iNodeDof));

                for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
                {
                    int dofID = nodePtr->GetDof(dof, iDof);
                    std::unordered_map<int, int>::const_iterator got = mMappingGlobal2LocalDof.find(dofID);
                    if (got == mMappingGlobal2LocalDof.end()) // not found => insert a new one
                    {
                        mMappingGlobal2LocalDof.insert(std::make_pair(dofID, localDofIndex));
                        //                        std::cout << "(dofID, localDofIndex): " << dofID << " " <<
                        //                        localDofIndex << std::endl;
                        localDofIndex++;
                        mNumMasterDofs++;
                    }
                }
            }
        }
    }

    mNumDofs = mMappingGlobal2LocalDof.size();
    mDofMappingComputed = true;
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillMappingGlobalLocalSlaveNodes()
{
    mNumSlaveNodes = 0;
    int numElements = 0;
    for (auto& it : mElementsSlave)
    {
        auto activeDofs = it.first->GetInterpolationType()->GetActiveDofs();
        if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
            throw MechanicsException(__PRETTY_FUNCTION__, "Contact element is only implemented for displacements.");

        auto dof = Node::eDof::DISPLACEMENTS;

        const InterpolationType* interpolationType = it.first->GetInterpolationType();
        const InterpolationBase& interpolationDof = interpolationType->Get(dof);
        const int numNodes = interpolationDof.GetNumNodes();

        numElements++;

        for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
        {
            const NodeBase* nodePtr = it.first->GetNode(interpolationDof.GetNodeIndex(iNodeDof));
            int globalNodeId = this->mStructure->NodeGetId(nodePtr);

            std::unordered_map<int, int>::const_iterator got = mMappingGlobal2LocalSlaveNodes.find(globalNodeId);
            if (got == mMappingGlobal2LocalSlaveNodes.end()) // not found => insert a new one
            {
                mMappingGlobal2LocalSlaveNodes.insert(std::make_pair(globalNodeId, mNumSlaveNodes));
                mNumSlaveNodes++;
            }
        }
    }

    mSlaveNodesMappingComputed = true;
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillMappingGlobalLocalMasterNodes()
{
    mNumMasterNodes = 0;

    for (int i = 0; i < mElementsMaster.rows(); i++)
    {
        for (int j = 0; j < mElementsMaster.cols(); j++)
        {
            const ContinuumElementIGA<TDimMaster>* masterElement = mElementsMaster(i, j).first;
            auto activeDofs = masterElement->GetInterpolationType()->GetActiveDofs();
            if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
                throw MechanicsException(__PRETTY_FUNCTION__, "Contact element is only implemented for displacements.");

            auto dof = Node::eDof::DISPLACEMENTS;

            const InterpolationType* interpolationType = masterElement->GetInterpolationType();
            const InterpolationBase& interpolationDof = interpolationType->Get(dof);
            const int numNodes = interpolationDof.GetNumNodes();

            for (int iNodeMaster = 0; iNodeMaster < numNodes; ++iNodeMaster)
            {
                const NodeBase* nodePtr = masterElement->GetNode(interpolationDof.GetNodeIndex(iNodeMaster));
                int globalNodeId = this->mStructure->NodeGetId(nodePtr);

                std::unordered_map<int, int>::const_iterator got = mMappingGlobal2LocalMasterNodes.find(globalNodeId);
                if (got == mMappingGlobal2LocalMasterNodes.end()) // not found => insert a new one
                {
                    mMappingGlobal2LocalMasterNodes.insert(std::make_pair(globalNodeId, mNumMasterNodes));
                    // std::cout << "(dofID, localDofIndex): " << dofID << " " << localDofIndex << std::endl;
                    mNumMasterNodes++;
                }
            }
        }
    }

    mMasterNodesMappingComputed = true;
}

namespace NuTo
{

template <>
const ContinuumContactElement<1, 1>& NuTo::ContinuumContactElement<1, 1>::AsContinuumContactElement11() const
{
    return *this;
}

template <>
ContinuumContactElement<1, 1>& NuTo::ContinuumContactElement<1, 1>::AsContinuumContactElement11()
{
    return *this;
}

template <>
const ContinuumContactElement<2, 1>& NuTo::ContinuumContactElement<2, 1>::AsContinuumContactElement21() const
{
    return *this;
}

template <>
ContinuumContactElement<2, 1>& NuTo::ContinuumContactElement<2, 1>::AsContinuumContactElement21()
{
    return *this;
}

template <>
const ContinuumContactElement<2, 2>& NuTo::ContinuumContactElement<2, 2>::AsContinuumContactElement22() const
{
    return *this;
}

template <>
ContinuumContactElement<2, 2>& NuTo::ContinuumContactElement<2, 2>::AsContinuumContactElement22()
{
    return *this;
}


template <>
ContinuumContactElement<3, 2>& NuTo::ContinuumContactElement<3, 2>::AsContinuumContactElement32()
{
    return *this;
}

} // end namespace NuTo

template class NuTo::ContinuumContactElement<3, 2>; // FEM <-> IGA L
template class NuTo::ContinuumContactElement<2, 2>; // FEM/IGA/IGA L <-> IGA/IGA L
template class NuTo::ContinuumContactElement<2, 1>; // FEM/IGA <-> IGA L
template class NuTo::ContinuumContactElement<1, 2>; // IGA L <-> IGA
template class NuTo::ContinuumContactElement<1, 1>; // IGA L <-> IGA L
