#include <boost/ptr_container/ptr_list.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <sstream>
#include <iostream>
#include <string>
#include "math/EigenSolverArpack.h"

#include "mechanics/structures/StructureBase.h"
#include "base/Timer.h"
#include "base/serializeStream/SerializeStreamIn.h"
#include "base/serializeStream/SerializeStreamOut.h"
#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseDirectSolverMKLPardiso.h"
#include "math/SparseDirectSolverPardiso.h"

#include "math/SparseMatrixCSRSymmetric.h"
#include "math/SparseMatrixCSRVector2General.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/groups/Group.h"
#include "mechanics/groups/GroupBase.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss13Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss16Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss3Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss6Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss12IpDetail.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NGauss9Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D4NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D4NGauss4Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NGauss1Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NGauss2x2x2Ip.h"
#include "mechanics/integrationtypes/IntegrationType3D8NLobatto.h"
#include "mechanics/integrationtypes/IntegrationType1D2NBoundaryGauss3Ip.h"
#include "mechanics/integrationtypes/IntegrationType0DBoundary.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/loads/LoadBase.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "base/Exception.h"
#include "mechanics/structures/StructureBaseEnum.h"
#include "mechanics/structures/StructureOutputBase.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/structures/StructureOutputBlockVector.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/structures/Assembler.h"
#include "mechanics/constraints/ConstraintCompanion.h"

#ifdef ENABLE_VISUALIZE
#include "visualize/UnstructuredGrid.h"
#endif // ENABLE_VISUALIZE
#include "visualize/ComponentName.h"

using namespace NuTo;

NuTo::StructureBase::StructureBase(int rDimension)
    : mAssembler(std::make_unique<Assembler>())
    , mShowTime(true)
    , mVerboseLevel(0)
{
    if (rDimension != 1 && rDimension != 2 && rDimension != 3)
    {
        throw Exception("[StructureBase::StructureBase] The dimension of a structure is either 1, 2 or 3.");
    }
    mDimension = rDimension;

    mNumTimeDerivatives = 0;

    mHaveTmpStaticData = false;
    mUpdateTmpStaticDataRequired = true;
    mToleranceStiffnessEntries = 0.;

#ifdef _OPENMP
    // then the environment variable is used
    mNumProcessors = 1;
#endif // _OPENMP
}

NuTo::StructureBase::~StructureBase()
{
}

int NuTo::StructureBase::GetDimension() const
{
    return mDimension;
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::StructureBase::Info() const
{
    mLogger << "dimension : " << mDimension << "\n";

    mLogger << "number of time derivatives : " << mNumTimeDerivatives << "\n";

    mLogger << "num dofs : " << GetNumTotalDofs() << "\n";
    mLogger << "num active dofs : " << GetNumTotalActiveDofs() << "\n";

    // print info for groups
    GroupInfo(mVerboseLevel);
}

//! @brief ... number of time derivatives for the nodes (0 : static, 1: velocities, 2: accelerations)
void NuTo::StructureBase::SetNumTimeDerivatives(int rNumTimeDerivatives)
{
    if (rNumTimeDerivatives < 0 || rNumTimeDerivatives > 2)
        throw NuTo::Exception(
                "[NuTo::StructureBase::SetNumTimeDerivatives] number of time derivatives is either 0, 1 or 2.");

    mNumTimeDerivatives = rNumTimeDerivatives;
}

//! @brief ... return number of time derivatives (0 : static, 1: velocities, 2: accelerations)
int NuTo::StructureBase::GetNumTimeDerivatives() const
{
    return mNumTimeDerivatives;
}


// store all elements of a group in a vector
void NuTo::StructureBase::GetElementsByGroup(Group<ElementBase>* rElementGroup, std::vector<ElementBase*>& rElements)
{
    Group<ElementBase>::iterator ElementIter = rElementGroup->begin();
    while (ElementIter != rElementGroup->end())
    {
        rElements.push_back(ElementIter->second);
        ++ElementIter;
    }
}

const NuTo::Constraint::Constraints& NuTo::StructureBase::Constraints() const
{
    return GetAssembler().GetConstraints();
}

NuTo::Constraint::Constraints& NuTo::StructureBase::Constraints()
{
    return GetAssembler().GetConstraints();
}

void NuTo::StructureBase::AddVisualizationComponent(int rElementGroup, eVisualizeWhat visualizeComponent)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Element group does not exist.");

    // create a new visualization list for an element group or add components to an already existing list
    if (mGroupVisualizeComponentsMap.find(rElementGroup) == mGroupVisualizeComponentsMap.end())
        mGroupVisualizationType.emplace(rElementGroup, eVisualizationType::VORONOI_CELL);

    mGroupVisualizeComponentsMap[rElementGroup].push_back(visualizeComponent);

    if (mVerboseLevel > 5)
    {
        for (auto const& iPair : mGroupVisualizeComponentsMap)
        {
            std::cout << "element group: \t" << iPair.first << std::endl;
            for (auto component : iPair.second)
            {
                std::cout << "components: \t " << GetComponentName(component) << std::endl;
            }
        }
    }
#endif // ENABLE_VISUALIZE
}


void NuTo::StructureBase::AddVisualizationComponent(int rElementGroup, const std::string& visualizeComponent)
{
    AddVisualizationComponent(rElementGroup, GetComponentEnum(visualizeComponent));
}

void NuTo::StructureBase::SetVisualizationType(const int rElementGroup, const eVisualizationType rVisualizationType)
{
    // check if the element group exists
    if (mGroupMap.find(rElementGroup) == mGroupMap.end())
        throw Exception(__PRETTY_FUNCTION__, "Element group does not exist.");

    // check if the element group exists
    if (mGroupVisualizationType.find(rElementGroup) == mGroupVisualizationType.end())
        throw Exception(__PRETTY_FUNCTION__,
                        "Please add a visualization component first before setting the visualization type.");

    mGroupVisualizationType.at(rElementGroup) = rVisualizationType;
}


void NuTo::StructureBase::ClearVisualizationComponents()
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    mGroupVisualizeComponentsMap.clear();

#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::ExportVtkDataFileNodes(const std::string& rResultFileName)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    for (auto const& it : mGroupVisualizeComponentsMap)
    {
        Visualize::UnstructuredGrid visualize;
        DefineVisualizeNodeData(visualize, it.second);
        NodeTotalAddToVisualize(visualize, it.second);

        visualize.ExportVtuDataFile(rResultFileName);
    }
#endif // ENABLE_VISUALIZE
}


void NuTo::StructureBase::ExportVtkDataFileElements(const std::string& rResultFileName, bool asBinary)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    for (auto const& it : mGroupVisualizeComponentsMap)
    {
        Visualize::UnstructuredGrid visualize;

        DefineVisualizeElementData(visualize, it.second);
        ElementGroupAddToVisualize(it.first, visualize, it.second);

        visualize.ExportVtuDataFile(rResultFileName, asBinary);
    }

#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::ElementGroupExportVtkDataFile(int rGroupIdent, const std::string& rResultFileName)
{
#ifdef ENABLE_VISUALIZE
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    Visualize::UnstructuredGrid visualizer;
    DefineVisualizeElementData(visualizer, mGroupVisualizeComponentsMap.at(rGroupIdent));
    ElementGroupAddToVisualize(rGroupIdent, visualizer, mGroupVisualizeComponentsMap.at(rGroupIdent));

    visualizer.ExportVtuDataFile(rResultFileName);

#endif // ENABLE_VISUALIZE
}

std::map<int, std::vector<NuTo::eVisualizeWhat>>& NuTo::StructureBase::GetGroupVisualizeComponentsMap(void)
{
    return mGroupVisualizeComponentsMap;
}

void NuTo::StructureBase::CalculateInitialValueRates(TimeIntegrationBase& rTimeIntegrationScheme)
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

void NuTo::StructureBase::DefineVisualizeElementData(Visualize::UnstructuredGrid& visualizer,
                                                     const std::vector<NuTo::eVisualizeWhat>& visualizeComponents) const
{
#ifdef ENABLE_VISUALIZE

    for (auto component : visualizeComponents)
    {
        switch (component)
        {
        case NuTo::eVisualizeWhat::TOTAL_INELASTIC_EQ_STRAIN:
        case NuTo::eVisualizeWhat::LOCAL_EQ_STRAIN:
        case NuTo::eVisualizeWhat::DAMAGE:
        case NuTo::eVisualizeWhat::HEAT_FLUX:
        case NuTo::eVisualizeWhat::SLIP:
        case NuTo::eVisualizeWhat::PRINCIPAL_ENGINEERING_STRESS:
        case NuTo::eVisualizeWhat::LATTICE_STRESS:
        case NuTo::eVisualizeWhat::LATTICE_STRAIN:
        case NuTo::eVisualizeWhat::LATTICE_PLASTIC_STRAIN:
        case NuTo::eVisualizeWhat::ELECTRIC_FIELD:
        case NuTo::eVisualizeWhat::ELECTRIC_DISPLACEMENT:
        case NuTo::eVisualizeWhat::BOND_STRESS:
        case NuTo::eVisualizeWhat::ENGINEERING_PLASTIC_STRAIN:
        case NuTo::eVisualizeWhat::ENGINEERING_STRAIN:
        case NuTo::eVisualizeWhat::ENGINEERING_STRESS:
        case NuTo::eVisualizeWhat::SHRINKAGE_STRAIN:
        case NuTo::eVisualizeWhat::THERMAL_STRAIN:
            visualizer.DefineCellData(GetComponentName(component));
            break;

        case NuTo::eVisualizeWhat::ELECTRIC_POTENTIAL:
        case NuTo::eVisualizeWhat::TEMPERATURE:
        case NuTo::eVisualizeWhat::RELATIVE_HUMIDITY:
        case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
        case NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION:
        case NuTo::eVisualizeWhat::CRACK_PHASE_FIELD:
        case NuTo::eVisualizeWhat::CRACK_PHASE_FIELD_VELOCITY:
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        case NuTo::eVisualizeWhat::VELOCITY:
        case NuTo::eVisualizeWhat::ACCELERATION:
            visualizer.DefinePointData(GetComponentName(component));
            break;

        case NuTo::eVisualizeWhat::PARTICLE_RADIUS:
        case NuTo::eVisualizeWhat::ROTATION:
        case NuTo::eVisualizeWhat::ANGULAR_VELOCITY:
        case NuTo::eVisualizeWhat::ANGULAR_ACCELERATION:
            // do nothing;
            break;

        default:
            throw Exception(__PRETTY_FUNCTION__, "undefined visualize component " + GetComponentName(component));
        }
    }
#endif // ENABLE_VISUALIZE
}

void NuTo::StructureBase::DefineVisualizeNodeData(Visualize::UnstructuredGrid& visualizer,
                                                  const std::vector<NuTo::eVisualizeWhat>& visualizeComponents) const
{
#ifdef ENABLE_VISUALIZE

    for (auto component : visualizeComponents)
    {
        switch (component)
        {
        case NuTo::eVisualizeWhat::DISPLACEMENTS:
        case NuTo::eVisualizeWhat::ROTATION:
        case NuTo::eVisualizeWhat::VELOCITY:
            // case NuTo::eVisualizeWhat::ACCELERATION:
            // case NuTo::eVisualizeWhat::ANGULAR_VELOCITY:
            // case NuTo::eVisualizeWhat::ANGULAR_ACCELERATION:
            // case NuTo::eVisualizeWhat::PARTICLE_RADIUS:
            // case NuTo::eVisualizeWhat::TEMPERATURE:
            // case NuTo::eVisualizeWhat::NONLOCAL_EQ_STRAIN:
            // case NuTo::eVisualizeWhat::RELATIVE_HUMIDITY:
            // case NuTo::eVisualizeWhat::WATER_VOLUME_FRACTION:
            visualizer.DefinePointData(GetComponentName(component));

        default:
            // do nothing for integration point data in the visualization list. However, the visualization of new dofs
            // needs to be added here!
            break;
        }
    }
#endif // ENABLE_VISUALIZE
}


NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian(eStructureOutput rOutput)
{
    Timer timer(std::string(__FUNCTION__) + ": " + StructureOutputToString(rOutput), GetShowTime(), GetLogger());
    NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    std::set<eStructureOutput> supportedTypes({eStructureOutput::HESSIAN0, eStructureOutput::HESSIAN1,
                                               eStructureOutput::HESSIAN2, eStructureOutput::HESSIAN2_LUMPED});
    if (supportedTypes.find(rOutput) == supportedTypes.end())
        throw Exception(__PRETTY_FUNCTION__,
                        StructureOutputToString(rOutput) + " is not a matrix type or is not supported right now.");

    StructureOutputBlockMatrix hessian(GetDofStatus(), true);

    std::map<eStructureOutput, StructureOutputBase*> evaluateMap;
    evaluateMap[rOutput] = &hessian;

    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);

    Evaluate(input, evaluateMap);

    return hessian;
}

NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian0()
{
    return BuildGlobalHessian(eStructureOutput::HESSIAN0);
}

NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian1()
{
    return BuildGlobalHessian(eStructureOutput::HESSIAN1);
}

NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian2()
{
    return BuildGlobalHessian(eStructureOutput::HESSIAN2);
}

NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian2Lumped()
{
    return BuildGlobalHessian(eStructureOutput::HESSIAN2_LUMPED);
}


NuTo::StructureOutputBlockVector NuTo::StructureBase::BuildGlobalInternalGradient()
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    StructureOutputBlockVector internalGradient(GetDofStatus(), true);

    std::map<eStructureOutput, StructureOutputBase*> evaluateMap;
    evaluateMap[eStructureOutput::INTERNAL_GRADIENT] = &internalGradient;

    ConstitutiveInputMap input;
    input[Constitutive::eInput::CALCULATE_STATIC_DATA] =
            std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::EULER_BACKWARD);

    Evaluate(input, evaluateMap);

    return internalGradient;
}

NuTo::StructureOutputBlockMatrix NuTo::StructureBase::BuildGlobalHessian0_CDF(double rDelta)
{
    Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    bool showTime = GetShowTime();
    SetShowTime(false);


    /*
     * TT:
     * Warning: this code is ugly, repetitive and error-prone.
     * This is mainly due to the J/K  and JJ, JK, KJ, KK stuff.
     * Please be very careful and test properly.
     */

    StructureOutputBlockMatrix hessian0_CDF(GetDofStatus(), true);
    auto internalGradient0 = BuildGlobalInternalGradient();
    auto dofValues = NodeExtractDofValues(0);


    for (auto dofCol : GetDofStatus().GetActiveDofTypes())
    {
        // modify active dof values --> entries JJ (R.J) and KJ (R.K)
        auto& columnDofValues = dofValues.J[dofCol];

        for (int iCol = 0; iCol < columnDofValues.rows(); ++iCol)
        {
            columnDofValues[iCol] += rDelta;
            NodeMergeDofValues(0, dofValues);
            auto internalGradient1 = BuildGlobalInternalGradient();

            StructureOutputBlockVector column = (internalGradient1 - internalGradient0) * (1 / rDelta);
            for (auto dofRow : GetDofStatus().GetActiveDofTypes())
            {
                // set JJ entries
                auto& colJ = column.J[dofRow];
                auto& matrixJJ = hessian0_CDF.JJ(dofRow, dofCol);
                for (int i = 0; i < colJ.rows(); ++i)
                    matrixJJ.AddValue(i, iCol, colJ(i, 0));

                // set KJ entries
                auto& colK = column.K[dofRow];
                auto& matrixKJ = hessian0_CDF.KJ(dofRow, dofCol);
                for (int i = 0; i < colK.rows(); ++i)
                    matrixKJ.AddValue(i, iCol, colK(i, 0));
            }
            columnDofValues[iCol] -= rDelta;
        }


        // modify dependent dof values --> entries JK (R.J) and KK (R.K)
        auto& rowDofValues = dofValues.K[dofCol];

        for (int iCol = 0; iCol < rowDofValues.rows(); ++iCol)
        {
            rowDofValues[iCol] += rDelta;
            NodeMergeDofValues(0, dofValues);
            auto internalGradient1 = BuildGlobalInternalGradient();

            StructureOutputBlockVector column = (internalGradient1 - internalGradient0) * (1 / rDelta);
            for (auto dofRow : GetDofStatus().GetActiveDofTypes())
            {
                // set JK entries
                auto& colJ = column.J[dofRow];
                auto& matrixJK = hessian0_CDF.JK(dofRow, dofCol);
                for (int i = 0; i < colJ.rows(); ++i)
                    matrixJK.AddValue(i, iCol, colJ(i, 0));

                // set KJ entries
                auto& colK = column.K[dofRow];
                auto& matrixKK = hessian0_CDF.KK(dofRow, dofCol);
                for (int i = 0; i < colK.rows(); ++i)
                    matrixKK.AddValue(i, iCol, colK(i, 0));
            }
            rowDofValues[iCol] -= rDelta;
        }
    }
    NodeMergeDofValues(0, dofValues);


    SetShowTime(showTime);

    return hessian0_CDF;
}

bool NuTo::StructureBase::CheckHessian0(double rDelta, double rRelativeTolerance, bool rPrintWrongMatrices)
{
    bool isHessianCorrect = true;

    NodeBuildGlobalDofs(__FUNCTION__);
    bool hasInteractingConstraints = GetDofStatus().HasInteractingConstraints();
    DofStatusSetHasInteractingConstraints(
            true); // this ensures the full assembly of KJ and KK, which could be skipped if CMat.Entries = 0

    auto hessian0 = BuildGlobalHessian0();
    auto hessian0_CDF = BuildGlobalHessian0_CDF(rDelta);

    isHessianCorrect = isHessianCorrect &&
                       CheckHessian0_Submatrix(hessian0.JJ, hessian0_CDF.JJ, rRelativeTolerance, rPrintWrongMatrices);
    isHessianCorrect = isHessianCorrect &&
                       CheckHessian0_Submatrix(hessian0.JK, hessian0_CDF.JK, rRelativeTolerance, rPrintWrongMatrices);
    isHessianCorrect = isHessianCorrect &&
                       CheckHessian0_Submatrix(hessian0.KJ, hessian0_CDF.KJ, rRelativeTolerance, rPrintWrongMatrices);
    isHessianCorrect = isHessianCorrect &&
                       CheckHessian0_Submatrix(hessian0.KK, hessian0_CDF.KK, rRelativeTolerance, rPrintWrongMatrices);

    DofStatusSetHasInteractingConstraints(hasInteractingConstraints);

    return isHessianCorrect;
}

bool NuTo::StructureBase::CheckHessian0_Submatrix(const BlockSparseMatrix& rHessian0, BlockSparseMatrix& rHessian0_CDF,
                                                  double rRelativeTolerance, bool rPrintWrongMatrices)
{
    int row, col;
    assert(rHessian0.GetNumRows() == rHessian0_CDF.GetNumRows());
    assert(rHessian0.GetNumColumns() == rHessian0_CDF.GetNumColumns());

    if (rHessian0.GetNumRows() == 0 or rHessian0.GetNumColumns() == 0)
        return true;

    bool isSubmatrixCorrect = true;
    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, " ", "\n", "|", " |");
    for (auto dofRow : GetDofStatus().GetActiveDofTypes())
    {
        for (auto dofCol : GetDofStatus().GetActiveDofTypes())
        {
            Eigen::MatrixXd hessian0_CDF_Full = rHessian0_CDF(dofRow, dofCol).ConvertToFullMatrix();

            double scaling = 1. / rHessian0_CDF(dofRow, dofCol).AbsMax();

            auto& diff = rHessian0_CDF(dofRow, dofCol);
            diff.AddScal(rHessian0(dofRow, dofCol), -1.);

            diff *= scaling;
            double error = diff.AbsMax(row, col);
            if (error > rRelativeTolerance)
            {
                GetLogger() << "[" << __FUNCTION__ << "] max error in (" << Node::DofToString(dofRow) << ","
                            << Node::DofToString(dofCol) << ") " << error << " at entry (" << row << "," << col
                            << ")\n";
                GetLogger() << "hessian0(" << row << "," << col
                            << ") = " << rHessian0(dofRow, dofCol).ConvertToFullMatrix()(row, col) << "\n";
                GetLogger() << "hessian0_CDF(" << row << "," << col << ") = " << hessian0_CDF_Full(row, col) << "\n";
                isSubmatrixCorrect = false;
                if (rPrintWrongMatrices)
                {
                    Eigen::MatrixXd diffPrint = diff.ConvertToFullMatrix();
                    GetLogger() << "####### relative difference\n" << diffPrint.format(fmt) << "\n";
                    GetLogger() << "####### hessian0\n"
                                << rHessian0(dofRow, dofCol).ConvertToFullMatrix().format(fmt) << "\n";
                    GetLogger() << "####### hessian0_CDF\n" << hessian0_CDF_Full.format(fmt) << "\n";
                }
            }
        }
    }
    return isSubmatrixCorrect;
}

void NuTo::StructureBase::SolveGlobalSystemStaticElastic()
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());

    if (GetNumTimeDerivatives() > 0)
        throw NuTo::Exception(__PRETTY_FUNCTION__, "Only use this method for a system with 0 time derivatives.");

    NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    StructureOutputBlockVector deltaDof_dt0(GetDofStatus(), true);
    deltaDof_dt0.J.SetZero();
    deltaDof_dt0.K = GetAssembler().GetConstraintRhs();

    auto hessian0 = BuildGlobalHessian0();

    auto residual = hessian0 * deltaDof_dt0 - BuildGlobalExternalLoadVector() + BuildGlobalInternalGradient();

    hessian0.ApplyCMatrix(GetAssembler().GetConstraintMatrix());
    residual.ApplyCMatrix(GetAssembler().GetConstraintMatrix());

    // reuse deltaDof_dt0
    deltaDof_dt0.J = SolveBlockSystem(hessian0.JJ, residual.J);

    deltaDof_dt0.K = NodeCalculateDependentDofValues(deltaDof_dt0.J);
    NodeMergeDofValues(0, deltaDof_dt0);
}

void NuTo::StructureBase::ConstraintLinearEquationNodeToElementCreate(int rNode, int rElementGroup,
                                                                      NuTo::Node::eDof rDofType,
                                                                      const double rTolerance,
                                                                      Eigen::Vector3d rNodeCoordOffset)
{
    const int dim = GetDimension();

    Eigen::VectorXd queryNodeCoords = NodeGetNodePtr(rNode)->Get(Node::eDof::COORDINATES);
    queryNodeCoords = queryNodeCoords + rNodeCoordOffset.head(dim);


    std::vector<int> elementGroupIds = GroupGetMemberIds(rElementGroup);

    ElementBase* elementPtr = nullptr;
    Eigen::VectorXd elementNaturalNodeCoords;
    bool nodeInElement = false;
    for (auto const& eleId : elementGroupIds)
    {
        elementPtr = ElementGetElementPtr(eleId);

        // Coordinate interpolation must be linear so the shape function derivatives are constant!
        assert(elementPtr->GetInterpolationType().Get(Node::eDof::COORDINATES).GetTypeOrder() ==
               Interpolation::eTypeOrder::EQUIDISTANT1);
        const Eigen::MatrixXd& derivativeShapeFunctionsGeometryNatural =
                elementPtr->GetInterpolationType()
                        .Get(Node::eDof::COORDINATES)
                        .DerivativeShapeFunctionsNatural(
                                Eigen::VectorXd::Zero(dim)); // just as _some_ point, as said, constant

        // real coordinates of every node in rElement
        Eigen::VectorXd elementNodeCoords = elementPtr->ExtractNodeValues(NuTo::Node::eDof::COORDINATES);

        switch (mDimension)
        {
        case 2:
        {
            Eigen::Matrix2d invJacobian =
                    dynamic_cast<ContinuumElement<2>*>(elementPtr)
                            ->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, elementNodeCoords)
                            .inverse();

            elementNaturalNodeCoords = invJacobian * (queryNodeCoords - elementNodeCoords.head(2));
        }
        break;
        case 3:
        {
            Eigen::Matrix3d invJacobian =
                    dynamic_cast<ContinuumElement<3>*>(elementPtr)
                            ->CalculateJacobian(derivativeShapeFunctionsGeometryNatural, elementNodeCoords)
                            .inverse();

            elementNaturalNodeCoords = invJacobian * (queryNodeCoords - elementNodeCoords.head(3));
        }
        break;

        default:
            throw NuTo::Exception(std::string(__PRETTY_FUNCTION__) + ": \t Only implemented for 2D and 3D");
        }


        if ((elementNaturalNodeCoords.array() > -rTolerance).all() and
            elementNaturalNodeCoords.sum() <= 1. + rTolerance)
        {
            nodeInElement = true;
            break;
        }
    }

    if (not nodeInElement)
    {
        GetLogger() << "Natural node coordinates: \n" << elementNaturalNodeCoords << "\n";
        throw Exception(__PRETTY_FUNCTION__, "Node is not inside any element.");
    }

    auto shapeFunctions =
            elementPtr->GetInterpolationType().Get(Node::eDof::DISPLACEMENTS).ShapeFunctions(elementNaturalNodeCoords);

    std::vector<Constraint::Equation> equations(dim); // default construction of Equation with rhs = Constant = 0
    for (int iDim = 0; iDim < dim; ++iDim)
    {
        equations[iDim].AddTerm(Constraint::Term(*NodeGetNodePtr(rNode), iDim, 1.));
    }


    for (int iNode = 0; iNode < shapeFunctions.rows(); ++iNode)
    {
        int localNodeId = elementPtr->GetInterpolationType().Get(Node::eDof::DISPLACEMENTS).GetNodeIndex(iNode);
        auto globalNode = elementPtr->GetNode(localNodeId, Node::eDof::DISPLACEMENTS);
        //        std::cout << "globalNodeId \t" << globalNodeId << std::endl;
        double coefficient = -shapeFunctions(iNode, 0);

        for (int iDim = 0; iDim < dim; ++iDim)
            equations[iDim].AddTerm(Constraint::Term(*globalNode, iDim, coefficient));
    }
    Constraints().Add(Node::eDof::DISPLACEMENTS, equations);
}

void NuTo::StructureBase::Contact(const std::vector<int>& rElementGroups)
{
}


NuTo::BlockFullVector<double> NuTo::StructureBase::SolveBlockSystem(const BlockSparseMatrix& rMatrix,
                                                                    const BlockFullVector<double>& rVector) const
{
    Eigen::VectorXd resultForSolver;
    std::unique_ptr<NuTo::SparseMatrixCSR<double>> matrixForSolver = rMatrix.ExportToCSR();

    //    try
    //    {
    //        EigenSolverArpack m;
    //        auto evs = m.GetSmallest(*matrixForSolver);
    //        auto evl = m.GetLargest(*matrixForSolver);
    //        std::cout << "EV smallest: " << evs.first << std::endl;
    //        std::cout << "EV largest:  " << evl.first << std::endl;
    //        std::cout << "EV Condition:   " << evl.first / evs.first << std::endl;
    //    }catch(...)
    //    {
    //        std::cout << "Error calculating EVs" << std::endl;
    //    }

    matrixForSolver->SetOneBasedIndexing();

// allocate solver
#if defined(HAVE_PARDISO) && defined(_OPENMP)
    NuTo::SparseDirectSolverPardiso mySolver(GetNumProcessors(), GetVerboseLevel()); // note: not the MKL version
#else
    NuTo::SparseDirectSolverMUMPS mySolver;
#endif

    mySolver.SetShowTime(GetShowTime());


    mySolver.Solve(*matrixForSolver, rVector.Export(), resultForSolver);

    return BlockFullVector<double>(-resultForSolver, GetDofStatus());
}


NuTo::StructureOutputBlockVector NuTo::StructureBase::BuildGlobalExternalLoadVector()
{
    NuTo::Timer timer(__FUNCTION__, GetShowTime(), GetLogger());
    NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    StructureOutputBlockVector externalLoad(GetDofStatus(), true);

    for (const auto& load : mLoadMap)
    {
        load.second->AddLoadToGlobalSubVectors(externalLoad);
    }
    return externalLoad;
}

//! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
//! values smaller than that one will not be added to the global matrix
void NuTo::StructureBase::SetToleranceStiffnessEntries(double rToleranceStiffnessEntries)
{
    mToleranceStiffnessEntries = rToleranceStiffnessEntries;
}

//! @brief absolute tolerance for entries of the global stiffness matrix (coefficientMatrix0)
//! values smaller than that one will not be added to the global matrix
double NuTo::StructureBase::GetToleranceStiffnessEntries() const
{
    return mToleranceStiffnessEntries;
}

std::set<NuTo::Node::eDof> NuTo::StructureBase::DofTypesGet() const
{
    return GetDofStatus().GetDofTypes();
}

void NuTo::StructureBase::DofTypeDeactivateAll()
{
    std::set<NuTo::Node::eDof> activeDofs = DofTypesGetActive();
    for (auto dof : activeDofs)
        DofTypeSetIsActive(dof, false);

    assert(DofTypesGetActive().empty());
}

void NuTo::StructureBase::DofTypeActivateAll()
{
    DofTypeSetIsActive(DofTypesGet());
}

void NuTo::StructureBase::DofTypeSetIsActive(Node::eDof rDofType, bool rIsActive)
{
    BOOST_FOREACH (auto interpolationTypePair, mInterpolationTypeMap)
    {
        auto& interpolationType = interpolationTypePair.second;
        if (interpolationType->IsDof(rDofType))
            interpolationType->SetIsActive(rIsActive, rDofType);
    }
    UpdateDofStatus();
}

void NuTo::StructureBase::DofTypeSetIsActive(const std::set<Node::eDof>& rActiveDofTypes)
{
    DofTypeDeactivateAll();
    for (auto activeDofType : rActiveDofTypes)
        DofTypeSetIsActive(activeDofType, true);
}

bool NuTo::StructureBase::DofTypeIsActive(Node::eDof rDofType) const
{
    const auto& activeDofTypes = DofTypesGetActive();
    return activeDofTypes.find(rDofType) != activeDofTypes.end();
}

void NuTo::StructureBase::DofTypeSetIsConstitutiveInput(Node::eDof rDofType, bool rIsConstitutiveInput)
{
    BOOST_FOREACH (auto interpolationTypePair, mInterpolationTypeMap)
    {
        auto& interpolationType = interpolationTypePair.second;
        if (interpolationType->IsDof(rDofType))
            interpolationType->SetIsConstitutiveInput(rIsConstitutiveInput, rDofType);
    }
}

void NuTo::StructureBase::DofTypeSetIsSymmetric(Node::eDof rDofType, bool rIsSymmetric)
{
    GetAssembler().mDofStatus.SetIsSymmetric(rDofType, rIsSymmetric);
}

bool NuTo::StructureBase::DofTypeIsSymmetric(Node::eDof rDofType) const
{
    return GetDofStatus().IsSymmetric(rDofType);
}

const NuTo::DofStatus& NuTo::StructureBase::GetDofStatus() const
{
    return GetAssembler().mDofStatus;
}

void NuTo::StructureBase::UpdateDofStatus()
{
    std::set<Node::eDof> dofTypes;
    std::set<Node::eDof> activeDofTypes;
    for (auto interpolationTypePair : mInterpolationTypeMap)
    {
        const std::set<Node::eDof>& dofs = interpolationTypePair.second->GetDofs();
        dofTypes.insert(dofs.begin(), dofs.end());

        const std::set<Node::eDof>& activeDofs = interpolationTypePair.second->GetActiveDofs();
        activeDofTypes.insert(activeDofs.begin(), activeDofs.end());
    }
    dofTypes.erase(Node::eDof::COORDINATES);
    activeDofTypes.erase(Node::eDof::COORDINATES);

    GetAssembler().mDofStatus.SetDofTypes(dofTypes);
    GetAssembler().mDofStatus.SetActiveDofTypes(activeDofTypes);

    GetAssembler().mDofStatus.SetHasInteractingConstraints(GetAssembler().GetConstraintMatrix().GetNumActiveEntires() !=
                                                           0);
}


void NuTo::StructureBase::DofStatusSetHasInteractingConstraints(bool rHasInteractingConstraints)
{
    GetAssembler().mDofStatus.SetHasInteractingConstraints(rHasInteractingConstraints);
}


int NuTo::StructureBase::GetNumTotalDofs() const
{
    return GetNumTotalActiveDofs() + GetNumTotalDependentDofs();
}

int NuTo::StructureBase::GetNumTotalActiveDofs() const
{
    int numTotalActiveDofs = 0;
    for (auto pair : GetDofStatus().GetNumActiveDofsMap())
        numTotalActiveDofs += pair.second;
    return numTotalActiveDofs;
}

int NuTo::StructureBase::GetNumTotalDependentDofs() const
{
    int numTotalActiveDofs = 0;
    for (auto pair : GetDofStatus().GetNumDependentDofsMap())
        numTotalActiveDofs += pair.second;
    return numTotalActiveDofs;
}

std::set<NuTo::Node::eDof> NuTo::StructureBase::DofTypesGetActive() const
{
    return GetDofStatus().GetActiveDofTypes();
}


int NuTo::StructureBase::GetNumDofs(Node::eDof rDofType) const
{
    return GetNumActiveDofs(rDofType) + GetNumDependentDofs(rDofType);
}

int NuTo::StructureBase::GetNumActiveDofs(Node::eDof rDofType) const
{
    auto it = GetDofStatus().GetNumActiveDofsMap().find(rDofType);
    if (it == GetDofStatus().GetNumActiveDofsMap().end())
        throw NuTo::Exception(std::string("[") + __PRETTY_FUNCTION__ + "] There are no " + Node::DofToString(rDofType) +
                              " dofs.");
    return it->second;
}

int NuTo::StructureBase::GetNumDependentDofs(Node::eDof rDofType) const
{
    auto it = GetDofStatus().GetNumDependentDofsMap().find(rDofType);
    if (it == GetDofStatus().GetNumDependentDofsMap().end())
        throw NuTo::Exception(std::string("[") + __PRETTY_FUNCTION__ + "] There are no " + Node::DofToString(rDofType) +
                              " dofs.");
    return it->second;
}

int NuTo::StructureBase::GetNumDofs(std::string rDofType) const
{
    return GetNumDofs(Node::DofToEnum(rDofType));
}

int NuTo::StructureBase::GetNumActiveDofs(std::string rDofType) const
{
    return GetNumActiveDofs(Node::DofToEnum(rDofType));
}

int NuTo::StructureBase::GetNumDependentDofs(std::string rDofType) const
{
    return GetNumDependentDofs(Node::DofToEnum(rDofType));
}

void NuTo::StructureBase::DofTypeSetIsActive(std::string rDofType, bool rIsActive)
{
    DofTypeSetIsActive(Node::DofToEnum(rDofType), rIsActive);
}

void NuTo::StructureBase::DofTypeSetIsConstitutiveInput(std::string rDofType, bool rIsConstitutiveInput)
{
    DofTypeSetIsConstitutiveInput(Node::DofToEnum(rDofType), rIsConstitutiveInput);
}

void NuTo::StructureBase::WriteRestartFile(std::string filename, double globalTime)
{
    NuTo::SerializeStreamOut binaryOut(filename, true);
    binaryOut << globalTime;
    binaryOut.Separator();
    binaryOut << *this;
}

double NuTo::StructureBase::ReadRestartFile(std::string filename)
{
    NuTo::SerializeStreamIn binaryIn(filename, true);
    double globalTime = -61.74;
    binaryIn >> globalTime;
    binaryIn.Separator();
    binaryIn >> *this;
    return globalTime;
}

#ifdef _OPENMP
//@brief determines the maximum independent sets and stores it at the structure
void NuTo::StructureBase::CalculateMaximumIndependentSets()
{
#define UNDONE 1
#define SELECTED 2
#define DELETED 3
    NuTo::Timer timer(__PRETTY_FUNCTION__, GetShowTime(), GetLogger());

    mMIS.clear();
    std::vector<ElementBase*> elementVector;
    GetElementsTotal(elementVector);

    // Build the connectivity graph
    // First get for all nodes all the elements
    std::map<const NodeBase*, std::vector<unsigned int>> elementsPerNode;
    for (unsigned int elementCount = 0; elementCount < elementVector.size(); elementCount++)
    {
        for (int nodeCount = 0; nodeCount < elementVector[elementCount]->GetNumInfluenceNodes(); nodeCount++)
        {
            elementsPerNode[elementVector[elementCount]->GetInfluenceNode(nodeCount)].push_back(elementCount);
        }
    }

    // Get the neighboring elements (always referring to the location in the vector elementVector)
    std::vector<std::vector<int>> NeighborElements(elementVector.size());
    for (auto& node : elementsPerNode)
    {
        for (unsigned int elementCount1 = 0; elementCount1 < node.second.size(); elementCount1++)
        {
            for (unsigned int elementCount2 = elementCount1 + 1; elementCount2 < node.second.size(); elementCount2++)
            {
                NeighborElements[node.second[elementCount1]].push_back(node.second[elementCount2]);
                NeighborElements[node.second[elementCount2]].push_back(node.second[elementCount1]);
            }
        }
    }

    // build the maximum independent sets
    std::vector<int> elementState(elementVector.size());
    for (int& entry : elementState)
        entry = UNDONE;

    unsigned int numDeleted = 0;
    unsigned int curMIS = 0;
    mMIS.resize(10);
    while (numDeleted < elementVector.size())
    {
        if (mMIS.size() <= curMIS)
            mMIS.resize(curMIS + 1);
        for (unsigned int countElement = 0; countElement < elementVector.size(); countElement++)
        {
            if (elementState[countElement] != UNDONE)
                continue;

            // add element to the set
            (mMIS[curMIS]).push_back(elementVector[countElement]);
            elementState[countElement] = DELETED;
            numDeleted++;

            // mark all the neighboring elements as selected, which prevents them to being added to this set
            for (int theNeighborCount : NeighborElements[countElement])
            {
                if (elementState[theNeighborCount] == UNDONE)
                    elementState[theNeighborCount] = SELECTED;
            }
        }
        // reset the selected elements to be undone
        for (unsigned int countElement = 0; countElement < elementVector.size(); countElement++)
        {
            if (elementState[countElement] == SELECTED)
                elementState[countElement] = UNDONE;
        }
        curMIS++;
    }
    mMIS.resize(curMIS);
}
#else
//@brief determines the maximum independent sets and stores it at the structure, do nothing for applications without
// openmp
void NuTo::StructureBase::CalculateMaximumIndependentSets()
{
}
#endif

//@brief set the number of processors for openmp parallelization
void NuTo::StructureBase::SetNumProcessors(int rNumProcessors)
{
#ifdef _OPENMP
    mNumProcessors = rNumProcessors;
#endif //_OPENMP
}
//@brief get the number of processors for openmp parallelization
int NuTo::StructureBase::GetNumProcessors() const
{
#ifdef _OPENMP
    return mNumProcessors;
#endif //_OPENMP
    return 1;
}

void NuTo::StructureBase::SetOMPNested(bool rNested)
{
#ifdef _OPENMP
    omp_set_nested(rNested);
#endif //_OPENMP
}

bool NuTo::StructureBase::InterpolationTypeIsConstitutiveInput(NuTo::Node::eDof rDofType)
{
    for (auto interpolation = mInterpolationTypeMap.begin(); interpolation != mInterpolationTypeMap.end();
         interpolation++)
    {
        auto interpolationType = interpolation->second;
        if (interpolationType->IsConstitutiveInput(rDofType))
        {
            return true;
        }
    }

    return false;
}

bool StructureBase::GetShowTime() const
{
    return mShowTime;
}

void StructureBase::SetShowTime(bool showTime)
{
    mShowTime = showTime;
}

unsigned short StructureBase::GetVerboseLevel() const
{
    return mVerboseLevel;
}

void StructureBase::SetVerboseLevel(unsigned short verboseLevel)
{
    mVerboseLevel = verboseLevel;
}
