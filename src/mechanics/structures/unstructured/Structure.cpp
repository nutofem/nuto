#include <boost/assign/ptr_map_inserter.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/spirit/include/classic_core.hpp>

#include "mechanics/structures/unstructured/Structure.h"

#include "base/Timer.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"

#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputDummy.h"
#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorInt.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/groups/GroupBase.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/timeIntegration/TimeIntegrationBase.h"

#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

#include "mechanics/mesh/MeshCompanion.h"
#include "mechanics/structures/Assembler.h"

NuTo::Structure::Structure(int rDimension)
    : StructureBase(rDimension)
{
}

NuTo::Structure::~Structure()
{
    mElementMap.clear();
}

void NuTo::Structure::Info() const
{
    StructureBase::Info();
    NodeInfo(mVerboseLevel);
    ElementInfo(mVerboseLevel);
}


void NuTo::Structure::Evaluate(const NuTo::ConstitutiveInputMap& rInput,
                               std::map<eStructureOutput, StructureOutputBase*>& rStructureOutput)
{
    std::string outputs = " ";
    for (auto it : rStructureOutput)
        outputs += StructureOutputToString(it.first) + " ";

    Timer timer(std::string(__FUNCTION__) + outputs, GetShowTime(), GetLogger());

    if (rStructureOutput.empty())
        return; // ! ---> may occur if matrices have been identified as constant

    NodeBuildGlobalDofs();

    // build global tmp static data
    if (this->mHaveTmpStaticData && this->mUpdateTmpStaticDataRequired)
        throw Exception(__PRETTY_FUNCTION__, "First update of tmp static data required.");

    for (auto iteratorOutput : rStructureOutput)
    {
        iteratorOutput.second->SetZero();
    }

#ifdef _OPENMP
    std::string exceptionMessage = "";
    if (mNumProcessors != 0)
    {
        omp_set_num_threads(mNumProcessors);
    }

    if (mMIS.size() == 0)
    {
        CalculateMaximumIndependentSets();
    }
    for (unsigned int misCounter = 0; misCounter < mMIS.size(); misCounter++)
    {
#pragma omp parallel shared(rStructureOutput) // firstprivate(elementOutputMap)
        {
#endif // _OPENMP
            // The allocation of the elementOutputMap is inside the openmp block
            // since the every thread needs a copy of the map.
            // This special case cannot (to my knowledge) be handled with the
            // omp firstprivate directive, since a copy of a shared_ptr is
            // not a deep copy of the underlying data - which makes perfectly sense.

            // BEWARE (!!!) Do not perform a SetZero on the rStructureOutput here
            // since it will remove the allocation of other MIS.

            std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> elementOutputMap;

            // allocate element outputs and resize the structure outputs
            for (auto iteratorOutput : rStructureOutput)
            {
                switch (iteratorOutput.first)
                {
                case NuTo::eStructureOutput::HESSIAN0:
                {
                    elementOutputMap[Element::eOutput::HESSIAN_0_TIME_DERIVATIVE] =
                            std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());
                    break;
                }
                case NuTo::eStructureOutput::HESSIAN1:
                {
                    elementOutputMap[Element::eOutput::HESSIAN_1_TIME_DERIVATIVE] =
                            std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());
                    break;
                }
                case NuTo::eStructureOutput::HESSIAN2:
                {
                    elementOutputMap[Element::eOutput::HESSIAN_2_TIME_DERIVATIVE] =
                            std::make_shared<ElementOutputBlockMatrixDouble>(GetDofStatus());
                    break;
                }
                case NuTo::eStructureOutput::HESSIAN2_LUMPED:
                {
                    elementOutputMap[Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE] =
                            std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());
                    break;
                }
                case NuTo::eStructureOutput::INTERNAL_GRADIENT:
                {
                    elementOutputMap[Element::eOutput::INTERNAL_GRADIENT] =
                            std::make_shared<ElementOutputBlockVectorDouble>(GetDofStatus());
                    break;
                }
                case NuTo::eStructureOutput::UPDATE_STATIC_DATA:
                {
                    elementOutputMap[Element::eOutput::UPDATE_STATIC_DATA] = std::make_shared<ElementOutputDummy>();
                    break;
                }
                default:
                {
                    throw NuTo::Exception(std::string("[") + __PRETTY_FUNCTION__ +
                                                   std::string("] Output request not implemented."));
                }
                }
            }
            // calculate element contribution
            elementOutputMap[Element::eOutput::GLOBAL_ROW_DOF] =
                    std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());
            elementOutputMap[Element::eOutput::GLOBAL_COLUMN_DOF] =
                    std::make_shared<ElementOutputBlockVectorInt>(GetDofStatus());
#ifdef _OPENMP
            for (auto elementIter = this->mMIS[misCounter].begin(); elementIter != this->mMIS[misCounter].end();
                 elementIter++)
            {
#pragma omp single nowait
                {
                    ElementBase* elementPtr = *elementIter;
                    // in OpenMP, exceptions may not leave the parallel region
                    try
                    {
                        elementPtr->Evaluate(rInput, elementOutputMap);
                    }
                    catch (std::exception& e)
                    {
                        exceptionMessage = e.what();
                    }

#else
    for (auto elementIter : this->mElementMap)
    {
        ElementBase* elementPtr = elementIter->second;
        elementPtr->Evaluate(rInput, elementOutputMap);
#endif

                    const auto& elementVectorGlobalDofsRow =
                            elementOutputMap.at(Element::eOutput::GLOBAL_ROW_DOF)->GetBlockFullVectorInt();
                    const auto& elementVectorGlobalDofsColumn =
                            elementOutputMap.at(Element::eOutput::GLOBAL_COLUMN_DOF)->GetBlockFullVectorInt();

                    for (auto& iteratorOutput : rStructureOutput)
                    {
                        StructureOutputBase* structureOutput = iteratorOutput.second;

                        switch (iteratorOutput.first)
                        {
                        case NuTo::eStructureOutput::HESSIAN0:
                        {
                            const auto& elementMatrix =
                                    elementOutputMap.at(Element::eOutput::HESSIAN_0_TIME_DERIVATIVE)
                                            ->GetBlockFullMatrixDouble();
                            structureOutput->AsStructureOutputBlockMatrix().AddElementMatrix(
                                    elementPtr, elementMatrix, elementVectorGlobalDofsRow,
                                    elementVectorGlobalDofsColumn, mToleranceStiffnessEntries,
                                    GetDofStatus().HasInteractingConstraints());
                            break;
                        }
                        case NuTo::eStructureOutput::HESSIAN1:
                        {
                            const auto& elementMatrix =
                                    elementOutputMap.at(Element::eOutput::HESSIAN_1_TIME_DERIVATIVE)
                                            ->GetBlockFullMatrixDouble();
                            structureOutput->AsStructureOutputBlockMatrix().AddElementMatrix(
                                    elementPtr, elementMatrix, elementVectorGlobalDofsRow,
                                    elementVectorGlobalDofsColumn, mToleranceStiffnessEntries,
                                    GetDofStatus().HasInteractingConstraints());
                            break;
                        }

                        case NuTo::eStructureOutput::HESSIAN2:
                        {
                            const auto& elementMatrix =
                                    elementOutputMap.at(Element::eOutput::HESSIAN_2_TIME_DERIVATIVE)
                                            ->GetBlockFullMatrixDouble();
                            structureOutput->AsStructureOutputBlockMatrix().AddElementMatrix(
                                    elementPtr, elementMatrix, elementVectorGlobalDofsRow,
                                    elementVectorGlobalDofsColumn, mToleranceStiffnessEntries,
                                    true); // always calculate the KJ and KK
                            // since its most likely only needed once,
                            // and causes troubles in the test files.
                            break;
                        }

                        case NuTo::eStructureOutput::HESSIAN2_LUMPED:
                        {
                            const auto& elementVector =
                                    elementOutputMap.at(Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE)
                                            ->GetBlockFullVectorDouble();

                            structureOutput->AsStructureOutputBlockMatrix().AddElementVectorDiagonal(
                                    elementVector, elementVectorGlobalDofsRow, mToleranceStiffnessEntries);
                            break;
                        }

                        case NuTo::eStructureOutput::INTERNAL_GRADIENT:
                        {
                            const auto& elementVector = elementOutputMap.at(Element::eOutput::INTERNAL_GRADIENT)
                                                                ->GetBlockFullVectorDouble();

                            structureOutput->AsStructureOutputBlockVector().AddElementVector(
                                    elementVector, elementVectorGlobalDofsRow);
                            break;
                        }

                        case NuTo::eStructureOutput::UPDATE_STATIC_DATA:
                            break;

                        default:
                        {
                            throw NuTo::Exception(__PRETTY_FUNCTION__,
                                                           StructureOutputToString(iteratorOutput.first) +
                                                                   " requested but not implemented.");
                        }
                        }
                    }

#ifdef _OPENMP
                }
            } // end loop over elements
        } // end parallel region
    } // end loop over independent sets

    if (exceptionMessage != "")
        throw Exception(exceptionMessage);
#else
    } // end loop over elements
#endif
}


void NuTo::Structure::CalculateInitialValueRates(NuTo::TimeIntegrationBase& rTimeIntegrationScheme)
{
    assert(mNumTimeDerivatives == 1 && "Using this function for 0 time derivatives seems to make no sense. More than "
                                       "one time derivative is not implemented so far!");


    constexpr const unsigned int maxIterations = 20;


    NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    rTimeIntegrationScheme.CalculateStaticAndTimeDependentExternalLoad();


    // declare necessary variables
    StructureOutputBlockMatrix Hessian_1(GetDofStatus(), true);

    StructureOutputBlockVector delta_dof_dt1(GetDofStatus(), true);

    StructureOutputBlockVector dof_dt1(GetDofStatus(), true);
    StructureOutputBlockVector residual(GetDofStatus(), true);
    StructureOutputBlockVector trialResidual(GetDofStatus(), true);
    StructureOutputBlockVector intForce(GetDofStatus(), true);
    StructureOutputBlockVector extForce(GetDofStatus(), true);


    // declare and fill output map for structure
    std::map<eStructureOutput, StructureOutputBase*> StructureOutputsTrial;
    StructureOutputsTrial[eStructureOutput::INTERNAL_GRADIENT] = &intForce;
    StructureOutputsTrial[eStructureOutput::HESSIAN1] = &Hessian_1;


    dof_dt1 = NodeExtractDofValues(1);

    ConstitutiveInputMap StructureInputs;
    StructureInputs[Constitutive::eInput::CALCULATE_INITIALIZE_VALUE_RATES] = nullptr;

    Evaluate(StructureInputs, StructureOutputsTrial);

    extForce = rTimeIntegrationScheme.CalculateCurrentExternalLoad(0);

    residual = intForce - extForce;


    unsigned int iteration = 0;


    while (residual.J.CalculateInfNorm() > rTimeIntegrationScheme.GetToleranceResidual())
    {
        ++iteration;
        if (iteration > maxIterations)
            throw Exception(__PRETTY_FUNCTION__, "No convergence while solving for initial value rates!");
        trialResidual = intForce - extForce;


        trialResidual.J -= Hessian_1.JJ * delta_dof_dt1.J;

        delta_dof_dt1.J = SolveBlockSystem(Hessian_1.JJ, residual.J);

        dof_dt1.J += delta_dof_dt1.J;

        NodeMergeDofValues(1, dof_dt1.J, dof_dt1.K);

        Evaluate(StructureInputs, StructureOutputsTrial);

        residual = intForce - extForce;
    }
}


std::vector<std::pair<int, int>> NuTo::Structure::ImportFromGmsh(const std::string& rFileName)
{
    return MeshCompanion::ImportFromGmsh(*this, rFileName, true);
}


std::vector<std::pair<int, int>> NuTo::Structure::ImportFromGmsh(const std::string& rFileName, bool useNewNumbers)
{
    return MeshCompanion::ImportFromGmsh(*this, rFileName, useNewNumbers);
}

std::vector<std::pair<int, int>> NuTo::Structure::ImportFromGmsh(const std::string &rFileName, bool useNewNumbers, std::map<int, int>& newNodes, std::map<int, int>& gmshNodes)
{
    return MeshCompanion::ImportFromGmsh(*this, rFileName, useNewNumbers, newNodes, gmshNodes);
}



void NuTo::Structure::CopyAndTranslate(Eigen::VectorXd& rOffset)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    std::map<NodeBase*, NodeBase*> old2NewNodePointer;
    std::map<ElementBase*, ElementBase*> old2NewElementPointer;
    CopyAndTranslate(rOffset, old2NewNodePointer, old2NewElementPointer);
}

void NuTo::Structure::CopyAndTranslate(Eigen::VectorXd& rOffset, std::map<NodeBase*, NodeBase*>& rOld2NewNodePointer,
                                       std::map<ElementBase*, ElementBase*>& rOld2NewElementPointer)
{
    if (rOffset.rows() != mDimension)
        throw Exception(__PRETTY_FUNCTION__, "offset has to have the same dimension as the structure.");
    if (rOffset.cols() != 1)
        throw Exception(__PRETTY_FUNCTION__, "offset has to have a single column.");

    std::vector<NodeBase*> nodeVector;
    GetNodesTotal(nodeVector);
    for (auto& node : nodeVector)
    {
        NodeBase* newNode = node->Clone();
        rOld2NewNodePointer[node] = newNode;

        // find unused integer id
        int id(mNodeMap.size());
        boost::ptr_map<int, NodeBase>::iterator it = mNodeMap.find(id);
        while (it != mNodeMap.end())
        {
            id++;
            it = mNodeMap.find(id);
        }

        // add node to map
        this->mNodeMap.insert(id, newNode);

        newNode->Set(Node::eDof::COORDINATES, node->Get(Node::eDof::COORDINATES) + rOffset);
    }
    // renumbering of dofs for global matrices required
    GetAssembler().SetNodeVectorChanged();

    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    for (auto oldElementPtr : elements)
    {
        int numNodes = oldElementPtr->GetNumNodes();
        std::vector<NodeBase*> nodeVector(numNodes);
        for (int countNode = 0; countNode < numNodes; countNode++)
        {
            nodeVector[countNode] = rOld2NewNodePointer[oldElementPtr->GetNode(countNode)];
        }
        // find interpolation type
        int interpolationTypeId = 0;
        const InterpolationType& interpolationTypeOld = oldElementPtr->GetInterpolationType();
        for (auto it = mInterpolationTypeMap.begin(); it != mInterpolationTypeMap.end(); it++)
            if ((it->second) == &interpolationTypeOld)
            {
                interpolationTypeId = it->first;
                break;
            }

        int newElementId = ElementCreate(interpolationTypeId, nodeVector);
        ElementBase* newElementPtr = ElementGetElementPtr(newElementId);
        rOld2NewElementPointer[oldElementPtr] = newElementPtr;

        // set integration type
        const IntegrationTypeBase& integrationType = oldElementPtr->GetIntegrationType();
        newElementPtr->SetIntegrationType(integrationType);

        // set section
        std::shared_ptr<const Section> section = oldElementPtr->GetSection();
        newElementPtr->SetSection(section);

        // set constitutive model
        ConstitutiveBase& constitutive = oldElementPtr->GetConstitutiveLaw(0);
        newElementPtr->SetConstitutiveLaw(constitutive);
    }
}

void NuTo::Structure::NuToSerializeSave(SerializeStreamOut& rStream)
{
    // be super carefull to symmetrically implement the same stuff to NuToSerializeLoad(...)

    // serialize nodes
    for (int i = 0; i < GetNumTimeDerivatives(); ++i)
    {
        auto nodalValues = NodeExtractDofValues(i);
        rStream << nodalValues.J;
        rStream.Separator();
        rStream << nodalValues.K;
        rStream.Separator();
    }

    // serialize element static data
    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    for (ElementBase* element : elements)
    {
        rStream << element->GetIPData();
        rStream.Separator();
    }
}

void NuTo::Structure::NuToSerializeLoad(SerializeStreamIn& rStream)
{
    // serialize nodes
    for (int i = 0; i < GetNumTimeDerivatives(); ++i)
    {
        auto nodalValues = NodeExtractDofValues(i);
        rStream >> nodalValues.J;
        rStream.Separator();
        rStream >> nodalValues.K;
        rStream.Separator();
        NodeMergeDofValues(i, nodalValues.J, nodalValues.K);
    }

    // serialize element static data
    std::vector<ElementBase*> elements;
    GetElementsTotal(elements);
    for (ElementBase* element : elements)
    {
        rStream >> element->GetIPData();
        rStream.Separator();
    }
}
