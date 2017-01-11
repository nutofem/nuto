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

template <int TDimSlave, int TDimMaster>
NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::ContinuumContactElement(const ContinuumElement<TDimSlave> *rSlaveElement, int rSurfaceId,
                                                                              Eigen::Matrix<std::pair<const ContinuumElementIGA<TDimMaster> *, int>, Eigen::Dynamic, Eigen::Dynamic> &rElementsMaster,
                                                                              const ConstitutiveBase *rConstitutiveContactLaw, int rContactAlgorithm)
    : ContinuumBoundaryElement<TDimSlave>(rSlaveElement, rSurfaceId), mElementsMaster(rElementsMaster), mNumDofs(0), mConstitutiveContactLaw(rConstitutiveContactLaw), mContactType(rContactAlgorithm)
{
    if(TDimMaster != ((rSurfaceId == -1) ? TDimSlave : TDimSlave - 1))
        throw MechanicsException(__PRETTY_FUNCTION__, "The dimension of master side interpolation is not correct.");

    mKnots.resize(TDimMaster);

    // initialize the knot values in x direction
    mKnots[0].resize(mElementsMaster.cols() + 1);
    int i = 0;
    for(; i < mElementsMaster.cols(); i++)
    {
        const ContinuumElementIGA<TDimMaster>* el = mElementsMaster(0,i).first;
        const Eigen::MatrixXd knots = el->GetKnots();
        mKnots[0](i) = knots(0,0);
    }
    mKnots[0](i) = mElementsMaster(0,i-1).first->GetKnots()(0,1);

    // initialize the knot values in y direction
    if(TDimMaster == 2)
    {
        int i = 0;
        for(; i < mElementsMaster.rows(); i++) mKnots[1](i) = mElementsMaster(i,0).first->GetKnots()(1,0);
        mKnots[1](i) = mElementsMaster(0,i-1).first->GetKnots()(1,1);
    }

    for(int dim = 0; dim < TDimMaster; dim++)
    {
        for(int i = 1; i < mKnots[dim].rows(); i++)
            if(mKnots[dim](i-1) > mKnots[dim](i))
                throw NuTo::MechanicsException(std::string(__PRETTY_FUNCTION__) + "Knots of inconsistent ordering.");
    }

    // TODO : check if ascending
    // check weather the elements are ordered
    for(int i = 1; i < mElementsMaster.rows(); i++)
    {
        for(int j = 1; j < mElementsMaster.cols(); j++)
        {
            mElementsMaster(i,j).first->GetKnots()(0,0);
        }
    }
}

template <int TDimSlave, int TDimMaster>
NuTo::Element::eElementType NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GetEnumType() const
{
    return Element::eElementType::CONTINUUMCONTACTELEMENT;
}


template<int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateGlobalRowDofs(BlockFullVector<int> &rGlobalRowDofs) const
{
    for (auto dof : this->mStructure->GetDofStatus().GetActiveDofTypes())
    {
        if (not (this->mInterpolationType->IsDof(dof)) && dof != Node::eDof::DISPLACEMENTS)
        {
            rGlobalRowDofs[dof].Resize(0);
            continue;
        }

        const InterpolationBase& interpolationType = this->mBaseElement->GetInterpolationType()->Get(dof);
        const int numNodes = interpolationType.GetNumNodes();

        FullVector<int, Eigen::Dynamic>& dofWiseGlobalRowDofs = rGlobalRowDofs[dof];

        int numSlaveDofs = interpolationType.GetNumDofs();
        dofWiseGlobalRowDofs.setZero(numSlaveDofs + mMappingGlobal2LocalDof.size());

        unsigned int numDofsPerType = this->mBaseElement->GetNode(interpolationType.GetNodeIndex(0))->GetNum(dof);

        for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
        {
            const NodeBase* nodePtr = this->mBaseElement->GetNode(interpolationType.GetNodeIndex(iNodeDof));

            for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
            {
                dofWiseGlobalRowDofs[numDofsPerType * iNodeDof + iDof] = nodePtr->GetDof(dof, iDof);
            }
        }

        // add master dofs
        for (auto& itMasterDofs: mMappingGlobal2LocalDof)
        {
            // itMasterDofs.second = local numbering
            // itMasterDofs.first = global numbering
            // both values are unique !!, see constructor
            dofWiseGlobalRowDofs[itMasterDofs.second] = itMasterDofs.first;
        }
    }
}

template<int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateGlobalColumnDofs(BlockFullVector<int> &rGlobalDofMapping) const
{
    if (this->GetNumNonlocalElements() == 0)
        CalculateGlobalRowDofs(rGlobalDofMapping);
    else
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented for nonlocal elements.");
}

template <int TDimSlave, int TDimMaster>
NuTo::eError NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::Evaluate(const ConstitutiveInputMap& rInput,
                                                                            std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput)
{
    EvaluateDataContinuumBoundary<TDimSlave> data;
    this->ExtractAllNecessaryDofValues(data);

    // the global dof numbering is not done in the constructor, because the dof numbering may not be completed at that time
    FillMappingGlobalLocal();

    auto constitutiveOutput = GetConstitutiveOutputMap(rElementOutput);
    auto constitutiveInput  = this->GetConstitutiveInputMap(constitutiveOutput);
    constitutiveInput.Merge(rInput);

    const InterpolationBase& interpolationType = this->mBaseElement->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
    int numNodesSlave = interpolationType.GetNumNodes();
    int numDofsSlave = interpolationType.GetNumDofs();
    int numAllDofsMaster = mMappingGlobal2LocalDof.size();

    data.mMortarGapMatrix.setZero(numDofsSlave + numAllDofsMaster, numNodesSlave);
    data.mMortarGapMatrixPenalty.setZero(numDofsSlave + numAllDofsMaster, numNodesSlave);
    data.mMortarGapVector.setZero(numNodesSlave, 1);
    data.mJacobianbyWeight.setZero(this->GetNumIntegrationPoints());

    for (int theIP = 0; theIP < this->GetNumIntegrationPoints(); theIP++)
    {
        this->CalculateNMatrixBMatrixDetJacobian(data, theIP);
        this->CalculateConstitutiveInputs(constitutiveInput, data);

        eError error = this->NuTo::ElementBase::EvaluateConstitutiveLaw<TDimSlave>(constitutiveInput, constitutiveOutput, theIP);

        if (error != eError::SUCCESSFUL)
            return error;

        for (auto it : rElementOutput)
        {
            switch (it.first)
            {
            case Element::eOutput::INTERNAL_GRADIENT:
            case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
                GapMatrixMortar(data, theIP);
                break;
            case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
            case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
            case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
            case Element::eOutput::GLOBAL_ROW_DOF:
            case Element::eOutput::GLOBAL_COLUMN_DOF:
            case Element::eOutput::UPDATE_STATIC_DATA:
            case Element::eOutput::UPDATE_TMP_STATIC_DATA:
                break;
            case Element::eOutput::IP_DATA:
                this->CalculateElementOutputIpData(it.second->GetIpData(), constitutiveOutput, theIP);
                break;
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "element output not implemented.");
            }
        }
    }

    CalculateElementOutputs(rElementOutput, data, constitutiveInput, constitutiveOutput);

    return eError::SUCCESSFUL;
}


template<int TDimSlave, int TDimMaster>
NuTo::ConstitutiveOutputMap NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GetConstitutiveOutputMap(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput) const
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
            this->FillConstitutiveOutputMapIpData(constitutiveOutput, it.second->GetIpData());
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
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillConstitutiveOutputMapHessian1(ConstitutiveOutputMap &rConstitutiveOutput,
                                                                                             BlockFullMatrix<double> &rHessian1) const
{
    (void) rConstitutiveOutput;
    auto activeDofs = this->mInterpolationType->GetActiveDofs();
    if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact is only implemented for displacements.");
    rHessian1(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS).Resize(mNumDofs, mNumDofs);
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillConstitutiveOutputMapHessian2Lumped(ConstitutiveOutputMap &rConstitutiveOutput,
                                                                                                   BlockFullVector<double> &rHessian2Lumped) const
{
    (void) rConstitutiveOutput;
    auto activeDofs = this->mInterpolationType->GetActiveDofs();
    if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact is only implemented for displacements.");
    rHessian2Lumped[Node::eDof::DISPLACEMENTS].Resize(mNumDofs);
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillConstitutiveOutputMapHessian2(ConstitutiveOutputMap &rConstitutiveOutput,
                                                                                                   BlockFullMatrix<double> &rHessian2) const
{
    (void) rConstitutiveOutput;
    auto activeDofs = this->mInterpolationType->GetActiveDofs();
    if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact is only implemented for displacements.");
    rHessian2(Node::eDof::DISPLACEMENTS,Node::eDof::DISPLACEMENTS).Resize(mNumDofs, mNumDofs);
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GapMatrixMortar(EvaluateDataContinuumBoundary<TDimSlave> &rData, int rTheIP)
{
    Eigen::VectorXd coordinatesIPSlave;
    Eigen::VectorXd ipCoordsNaturalSlave;
    GetGlobalIntegrationPointCoordinatesAndParameters(rTheIP, coordinatesIPSlave, ipCoordsNaturalSlave);

    // TODO mind the iga jacobian matrix (in case of iga layer- iga layer contact)
    auto jacobianSurface = this->mBaseElement->CalculateJacobianSurface(ipCoordsNaturalSlave, this->ExtractNodeValues(0, Node::eDof::COORDINATES), this->mSurfaceId);
    double jacobianbyWeight = jacobianSurface.norm() * this->GetIntegrationPointWeight(rTheIP);

    /*********************************************************************************/
    //  Projection of the rTheIP on the master element => \xi^s_{IP}, \xi^m_*, n^m_* //
    /*********************************************************************************/

    // **** Get the starting point for iteration **** //

    int numStartingPointsOneDir = 1;
    Eigen::VectorXd coords(numStartingPointsOneDir);
    coords << 0;

    int numStartingPointsElement = std::pow(numStartingPointsOneDir, TDimMaster);
    std::vector<Eigen::Matrix<double, TDimMaster, 1>> referenceCoordinates(numStartingPointsElement);

    int count = 0;
    for(int x = 0; x < numStartingPointsOneDir; x++)
    {
        referenceCoordinates[count](0,0) = coords(x);
        count++;
    }

    double minDistance = std::numeric_limits<double>::infinity();
    Eigen::VectorXd parameterMinMaster;
    Eigen::Vector2d indexMasterElement;
    for(int i = 0; i < mElementsMaster.rows(); i++)
    {
        for(int j = 0; j < mElementsMaster.cols(); j++)
        {
            const auto* elementPtr = mElementsMaster(i,j).first;
            int surfaceId = mElementsMaster(i,j).second;

            // ===> Get the position on the master curve/surface
            // ===> Compare and set to minimum if though
            const InterpolationBase& interpolationTypeDispMaster = elementPtr->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
            for(auto it : referenceCoordinates)
            {
                Eigen::VectorXd parameter = interpolationTypeDispMaster.CalculateNaturalSurfaceCoordinates(it, surfaceId, elementPtr->GetKnots());
                Eigen::VectorXd coordinatesMaster = elementPtr->InterpolateDofGlobalSurfaceDerivative(0, parameter, 0, 0);

                double distance = (coordinatesMaster - coordinatesIPSlave).norm();
                if(minDistance > distance)
                {
                    minDistance = distance;
                    parameterMinMaster = parameter;
                    indexMasterElement(0) = i;
                    indexMasterElement(1) = j;
                }
            }
        }
    }

    const ContinuumElementIGA<TDimMaster> *masterElement = mElementsMaster(indexMasterElement(0),indexMasterElement(1)).first;

    // **** Newton iteration ****//
    double tol = 1.e-10;
    double error = 1.;
    int maxNumIter = 100;
    int numIter = 0;
    Eigen::VectorXd r;
    Eigen::VectorXd coordinatesMaster;
    while(error > tol && numIter < maxNumIter)
    {
        // ==> function (dprime)
        coordinatesMaster = masterElement->InterpolateDofGlobalSurfaceDerivative(0, parameterMinMaster, 0, 0);
        r = coordinatesIPSlave - coordinatesMaster;
        Eigen::VectorXd dprime(TDimSlave - 1, 1);
        dprime.setZero(TDimSlave - 1);
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> prime = masterElement->InterpolateDofGlobalSurfaceDerivativeTotal(0, parameterMinMaster, 1);
        Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> primeprime = masterElement->InterpolateDofGlobalSurfaceDerivativeTotal(0, parameterMinMaster, 2);

        for(int j = 0; j < prime.cols(); j++)
            dprime(j) += r.dot(prime(0,j));

        // ==> derivative
        Eigen::MatrixXd dprimeprime(TDimSlave - 1, TDimSlave - 1);
        dprimeprime.setZero(TDimSlave - 1, TDimSlave - 1);

        for(int i = 0; i < prime.cols() ; i++)
        {
            for(int j = 0; j < prime.cols() ; j++)
            {
                dprimeprime(i,j) = r.dot(primeprime(i,j)) - prime(0,i).dot(prime(0,j));
            }
        }

        // iteration step
        Eigen::VectorXd increment = dprimeprime.colPivHouseholderQr().solve(-dprime);
        parameterMinMaster += increment;

        // New parameter value (check in which element)
        // for the FindSpan function degree = 0, since no multiple knots at the beginning and end
        if(TDimMaster == 1)
            indexMasterElement(1) = ShapeFunctionsIGA::FindSpan(parameterMinMaster(0), 0, mKnots[0]);
        else
        {
            indexMasterElement(0) = ShapeFunctionsIGA::FindSpan(parameterMinMaster(0), 0, mKnots[0]);
            indexMasterElement(1) = ShapeFunctionsIGA::FindSpan(parameterMinMaster(1), 0, mKnots[1]);
        }

        masterElement = mElementsMaster(indexMasterElement(0),indexMasterElement(1)).first;

        error = dprime.norm();
        numIter++;
    }

    if(numIter >= maxNumIter) std::cout << "!!!!!!ContinuumContactElement: Maximum number of Newton iterations exceeded!" << std::endl;


    /*************************/
    //  Build the gap matrix //
    /*************************/

    // normal vector
    Eigen::VectorXd normal = masterElement->InterpolateDofGlobalSurfaceNormal(parameterMinMaster);

    // Assemble the element gap matrix => \int_{\Gamma_e} F(ShapeFunctionsSlave(\xi^s), ShapeFunctionsMaster(\xi^*), n^*) d\Gamma

    const InterpolationBase& interpolationTypeCoordsSlave = this->mBaseElement->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
    Eigen::VectorXd  shapeFunsSlave = interpolationTypeCoordsSlave.CalculateShapeFunctions(ipCoordsNaturalSlave);

    const InterpolationBase& interpolationTypeCoordsMaster  = masterElement->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);
    Eigen::VectorXd  shapeFunsMaster = interpolationTypeCoordsMaster.CalculateShapeFunctions(parameterMinMaster); // master is always iga

    int numSlaveFunsNormal  = shapeFunsSlave.rows()*normal.rows();
    int numMasterFunsNormal = shapeFunsMaster.rows()*normal.rows();

    Eigen::VectorXd shapeFunsSlaveNormal(numSlaveFunsNormal);
    count = 0;
    for(int i = 0; i < shapeFunsSlave.rows(); i++)
    {
        shapeFunsSlaveNormal.block(count, 0, normal.rows(), 1) = shapeFunsSlave(i)*normal;
        count+=normal.rows();
    }

    Eigen::VectorXd shapeFunsMasterNormal(numMasterFunsNormal);
    count = 0;
    for(int i = 0; i < shapeFunsMaster.rows(); i++)
    {
        shapeFunsMasterNormal.block(count, 0, normal.rows(), 1) = (-1)*shapeFunsMaster(i)*normal;
        count+=normal.rows();
    }

    // matrix containing the derivatives of the slave and master side multiplied by the normal (R^S * normal \\  R^M * normal)
    Eigen::MatrixXd NContact;
    NContact.resize(numSlaveFunsNormal + numMasterFunsNormal, shapeFunsSlave.rows());

    NContact.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows())  =  shapeFunsSlaveNormal*shapeFunsSlave.transpose();
    NContact.block(numSlaveFunsNormal, 0, numMasterFunsNormal, shapeFunsSlave.rows()) = shapeFunsMasterNormal * shapeFunsSlave.transpose();
    NContact *= jacobianbyWeight;

    double gap = r.dot(normal);
//    if(gap < 0)
//        std::cout << gap << std::endl;

    double derivativeContactForce = mConstitutiveContactLaw->GetContactForceDerivative(gap);

    rData.mMortarGapMatrix.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows())        += NContact.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows());
    rData.mMortarGapMatrixPenalty.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows()) += NContact.block(0, 0, numSlaveFunsNormal, shapeFunsSlave.rows())*derivativeContactForce;

    const int numNodes = interpolationTypeCoordsMaster.GetNumNodes();
    unsigned int numDofsPerType = masterElement->GetNode(interpolationTypeCoordsMaster.GetNodeIndex(0))->GetNum(Node::eDof::DISPLACEMENTS);

    count = 0;
    for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
    {
        const NodeBase* nodePtr = masterElement->GetNode(interpolationTypeCoordsMaster.GetNodeIndex(iNodeDof));

        for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
        {
            int index = mMappingGlobal2LocalDof[nodePtr->GetDof(Node::eDof::DISPLACEMENTS, iDof)];
            rData.mMortarGapMatrix.row(index) += NContact.row(count + numSlaveFunsNormal);
            rData.mMortarGapMatrixPenalty.row(index) += (NContact.row(count + numSlaveFunsNormal)*derivativeContactForce);
            count++;
        }
    }

    // *** penalty vector or gap vector for the contact force ***
    if(mContactType == 0)
    {
        //rData.mMortarGapVector += NContact.transpose()*positions;
        rData.mMortarGapVector += shapeFunsSlave*gap*jacobianbyWeight;
    }
    else if(mContactType == 1)
    {
        rData.mMortarGapVector += shapeFunsSlave*(mConstitutiveContactLaw->GetContactForce(gap))*jacobianbyWeight;
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateElementOutputs(std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> &rElementOutput,
                                                                                   EvaluateDataContinuumBoundary<TDimSlave>                       &rData,
                                                                                   const ConstitutiveInputMap                                     &constitutiveInput,
                                                                                   const ConstitutiveOutputMap                                    &constitutiveOutput) const
{

    for (auto it : rElementOutput)
    {
        switch (it.first)
        {
        case Element::eOutput::INTERNAL_GRADIENT:
            CalculateElementOutputContactForce(it.second->GetBlockFullVectorDouble(), rData, constitutiveInput, constitutiveOutput);
            break;
        case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
            CalculateElementOutputContactForceDerivative(it.second->GetBlockFullMatrixDouble(), rData, constitutiveInput, constitutiveOutput);
            break;
        case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
        case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
        case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE: // already set to zero in GetConstitutiveOutputMap()...
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
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateElementOutputContactForce(BlockFullVector<double> &rInternalGradient,
                                                                                              EvaluateDataContinuumBoundary<TDimSlave> &rData,
                                                                                              const ConstitutiveInputMap &constitutiveInput,
                                                                                              const ConstitutiveOutputMap &constitutiveOutput) const
{
    for (auto dofRow : this->mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
        {
            if(mContactType == 0)
            {
                Eigen::VectorXd penaltyForce(rData.mMortarGapVector.rows());
                for(int i = 0; i < penaltyForce.rows(); i++)
                {
                    if (rData.mMortarGapVector(i) > 0.)
                        penaltyForce(i) = 0.;
                    else if (rData.mMortarGapVector(i) <= 0.)
                    {
                        penaltyForce(i) = mConstitutiveContactLaw->GetContactForce(rData.mMortarGapVector(i));
                    }
                }

                rInternalGradient[dofRow] = rData.mMortarGapMatrix*penaltyForce;
            }
            else if(mContactType == 1)
            {
                rInternalGradient[dofRow] = rData.mMortarGapMatrix*rData.mMortarGapVector;
            }
            break;
        }
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Element output CONTACT_FORCE for " + Node::DofToString(dofRow) + " not implemented.");
        }
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::CalculateElementOutputContactForceDerivative(BlockFullMatrix<double> &rGapMatrix,
                                                                                                        EvaluateDataContinuumBoundary<TDimSlave> &rData,
                                                                                                        const ConstitutiveInputMap &constitutiveInput,
                                                                                                        const ConstitutiveOutputMap &constitutiveOutput) const
{
    for (auto dofRow : this->mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : this->mInterpolationType->GetActiveDofs())
        {
            switch (Node::CombineDofs(dofRow, dofCol))
            {
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
            {
                rGapMatrix(dofRow, dofCol) = rData.mMortarGapMatrix*rData.mMortarGapMatrixPenalty.transpose();
//                if(rGapMatrix(dofRow, dofCol).norm() > 1.e-8)
//                    std::cout << "CCE: " << rGapMatrix(dofRow, dofCol).norm() << std::endl;
                break;
            }
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Element output CONTACT_FORCE_DERIVATIVE for " + Node::DofToString(dofRow) + " not implemented.");
            }
        }
    }
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GetGlobalIntegrationPointCoordinatesAndParameters(int rIpNum, Eigen::VectorXd &rCoordinatesIPSlave, Eigen::VectorXd &rParamsIPSlave) const
{
    Eigen::VectorXd naturalSurfaceIpCoordinates;
    switch (TDimMaster)
    {
        case 1:
        {
            double ipCoordinate;
            this->GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, ipCoordinate);
            naturalSurfaceIpCoordinates.resize(1);
            naturalSurfaceIpCoordinates(0) = ipCoordinate;
            break;
        }
        case 2:
        {
            double ipCoordinates[2];
            this->GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, ipCoordinates);
            naturalSurfaceIpCoordinates.resize(2);
            naturalSurfaceIpCoordinates(0) = ipCoordinates[0];
            naturalSurfaceIpCoordinates(1) = ipCoordinates[1];
            break;
        }
        default:
            break;
    }

    // ===> Get the position \xi^s_{IP}
    const InterpolationBase& interpolationTypeCoordsSlave = this->mBaseElement->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS);

    if (interpolationTypeCoordsSlave.GetTypeOrder() == Interpolation::eTypeOrder::SPLINE)
        rParamsIPSlave = interpolationTypeCoordsSlave.CalculateNaturalSurfaceCoordinates(naturalSurfaceIpCoordinates, this->mSurfaceId, this->mBaseElement->GetKnots()); // IGA
    else
        rParamsIPSlave = interpolationTypeCoordsSlave.CalculateNaturalSurfaceCoordinates(naturalSurfaceIpCoordinates, this->mSurfaceId); // FEM

    rCoordinatesIPSlave  = this->mBaseElement->InterpolateDofGlobalCurrentConfiguration(0, rParamsIPSlave, Node::eDof::COORDINATES, Node::eDof::DISPLACEMENTS);
}

template <int TDimSlave, int TDimMaster>
const Eigen::Vector3d NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::GetGlobalIntegrationPointCoordinates(int rIpNum) const
{
    Eigen::VectorXd naturalSurfaceIpCoordinates;
    switch (TDimMaster)
    {
        case 1:
        {
            double ipCoordinate;
            this->GetIntegrationType()->GetLocalIntegrationPointCoordinates1D(rIpNum, ipCoordinate);
            naturalSurfaceIpCoordinates.resize(1);
            naturalSurfaceIpCoordinates(0) = ipCoordinate;
            break;
        }
        case 2:
        {
            double ipCoordinates[2];
            this->GetIntegrationType()->GetLocalIntegrationPointCoordinates2D(rIpNum, ipCoordinates);
            naturalSurfaceIpCoordinates.resize(2);
            naturalSurfaceIpCoordinates(0) = ipCoordinates[0];
            naturalSurfaceIpCoordinates(1) = ipCoordinates[1];
            break;
        }
        default:
            break;
    }

    Eigen::VectorXd naturalIpCoordinates = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateNaturalSurfaceCoordinates(naturalSurfaceIpCoordinates, this->mSurfaceId);

    Eigen::MatrixXd matrixN = this->mInterpolationType->Get(Node::eDof::COORDINATES).CalculateMatrixN(naturalIpCoordinates);
    Eigen::VectorXd nodeCoordinates = this->ExtractNodeValues(0, Node::eDof::COORDINATES);

    Eigen::Vector3d globalIntegrationPointCoordinates = Eigen::Vector3d::Zero();
    globalIntegrationPointCoordinates.segment(0, this->GetLocalDimension()) = matrixN * nodeCoordinates;

    return globalIntegrationPointCoordinates;
}

template <int TDimSlave, int TDimMaster>
void NuTo::ContinuumContactElement<TDimSlave, TDimMaster>::FillMappingGlobalLocal()
{
    auto activeDofs = this->mInterpolationType->GetActiveDofs();
    if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "Contact element is only implemented for displacements.");

    int numSlaveDofs = this->mBaseElement->GetInterpolationType()->Get(Node::eDof::DISPLACEMENTS).GetNumDofs();

    int localDofIndex = numSlaveDofs;
    for(int i = 0; i < mElementsMaster.rows(); i++)
    {
        for(int j = 0; j < mElementsMaster.cols(); j++)
        {
            const ContinuumElementIGA<TDimMaster> *masterElement = mElementsMaster(i,j).first;
            auto activeDofs = masterElement->GetInterpolationType()->GetActiveDofs();
            if (activeDofs.size() > 1 && activeDofs.find(Node::eDof::DISPLACEMENTS) == activeDofs.end())
                throw MechanicsException(__PRETTY_FUNCTION__, "Contact element is only implemented for displacements.");

            auto dof = Node::eDof::DISPLACEMENTS;

            const InterpolationType *interpolationType = masterElement->GetInterpolationType();
            const InterpolationBase &interpolationDof = interpolationType->Get(dof);
            const int numNodes = interpolationDof.GetNumNodes();

            unsigned int numDofsPerType = masterElement->GetNode(interpolationDof.GetNodeIndex(0))->GetNum(dof);

            for (int iNodeDof = 0; iNodeDof < numNodes; ++iNodeDof)
            {
                const NodeBase* nodePtr = masterElement->GetNode(interpolationDof.GetNodeIndex(iNodeDof));

                for (unsigned iDof = 0; iDof < numDofsPerType; ++iDof)
                {
                    int dofID = nodePtr->GetDof(dof, iDof);
                    std::unordered_map<int,int>::const_iterator got = mMappingGlobal2LocalDof.find(dofID);
                    if (got == mMappingGlobal2LocalDof.end()) // not found => insert a new one
                    {
                        mMappingGlobal2LocalDof.insert( std::make_pair(dofID , localDofIndex) );
                        localDofIndex++;
                    }
                }
            }
        }
    }

    mNumDofs = numSlaveDofs + mMappingGlobal2LocalDof.size();
}

namespace NuTo
{

} // end namespace NuTo

template class NuTo::ContinuumContactElement<3,2>; // FEM <-> IGA L
template class NuTo::ContinuumContactElement<2,2>; // FEM/IGA/IGA L <-> IGA/IGA L
template class NuTo::ContinuumContactElement<2,1>; // FEM/IGA <-> IGA
template class NuTo::ContinuumContactElement<1,2>; // IGA L <-> IGA
template class NuTo::ContinuumContactElement<1,1>; // IGA L <-> IGA L

