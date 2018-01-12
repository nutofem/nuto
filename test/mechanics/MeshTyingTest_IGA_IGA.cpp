#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockVector.h"

#include <eigen3/Eigen/Core>

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverPardiso.h"

#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGaussNIp.h"

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#include "nuto/mechanics/IGA/NURBSCurve.h"
#include "nuto/mechanics/IGA/NURBSSurface.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#include "nuto/mechanics/groups/GroupEnum.h"

#include <boost/filesystem.hpp>

#include "nuto/mechanics/elements/ElementShapeFunctions.h"

#include "nuto/mechanics/timeIntegration/RungeKutta4.h"

#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockScalar.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include "nuto/mechanics/groups/GroupBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataBase.h"

#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include "nuto/mechanics/elements/ContinuumContactElement.h"

#define _USE_MATH_DEFINES
#include <cmath>

#ifdef ENABLE_VISUALIZE
#include "nuto/visualize/VisualizeEnum.h"
#endif

#define PRINTRESULT true

NuTo::NURBSSurface buildRect2DIGA(int rDegree, double x0, double y0, double Height, double Length,
                                  int &groupElements, int &groupNodes, int refinementX, int refinementY,
                                  NuTo::Structure *myStructure)
{
    /** Knots and control points **/
    int numElementsX = 1;
    int numElementsY = 1;

    Eigen::Vector2i degree(rDegree,rDegree);

    int numKnotsX = 2*(degree(0)+1) + numElementsX - 1;
    int numKnotsY = 2*(degree(1)+1) + numElementsY - 1;

    Eigen::VectorXd knotsX(numKnotsX);
    Eigen::VectorXd knotsY(numKnotsY);

    for (int i = 0; i <= degree(0); i++) knotsX(i) = 0.;
    for (int i = degree(0) + 1; i <= degree(0) + numElementsX - 1 ; i++) knotsX(i) = knotsX(i-1) + 1./numElementsX;
    for (int i = degree(0) + numElementsX ; i < numKnotsX; i++) knotsX(i) = 1;

    for (int i = 0; i <= degree(1); i++) knotsY(i) = 0.;
    for (int i = degree(1) + 1; i <= degree(1) + numElementsY - 1 ; i++) knotsY(i) = knotsY(i-1) + 1./numElementsY;
    for (int i = degree(1) + numElementsY ; i < numKnotsY; i++) knotsY(i) = 1;

    int numControlPointsX = (numKnotsX -1) - degree(0);
    int numControlPointsY = (numKnotsY -1) - degree(1);

    Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> controlPoints(numControlPointsY, numControlPointsX);

    double incrx = (Length/(numControlPointsX-1));
    double incry = (Height/(numControlPointsY-1));
    for(int i = 0; i < numControlPointsY; i++)
    {
        for(int j = 0; j < numControlPointsX; j++)
        {
            controlPoints(i,j) = Eigen::Vector2d(x0 + incrx*j, y0 + incry*i);
        }
    }

    Eigen::MatrixXd weights(numControlPointsY, numControlPointsX);
    weights.setOnes(numControlPointsY, numControlPointsX);

    NuTo::NURBSSurface surface(degree, knotsX, knotsY, controlPoints, weights);

    for(int i = 0; i < refinementX-1; i++) surface.DuplicateKnots(0);
    for(int i = 0; i < refinementY-1; i++) surface.DuplicateKnots(1);

    std::cout << "Control point:\t" << std::endl;
    for(int i = 0; i < surface.GetNumControlPoints(1); i++)
    {
        std::cout << std::endl;
        for(int j = 0; j < surface.GetNumControlPoints(0); j++)
            std::cout << "(" << surface.GetControlPoint(i,j).transpose() << ")\t";
    }

    std::cout << std::flush;

    std::set<NuTo::Node::eDof> setOfDOFS;
    setOfDOFS.insert(NuTo::Node::eDof::COORDINATES);
    setOfDOFS.insert(NuTo::Node::eDof::DISPLACEMENTS);

    surface.buildIGAStructure(*myStructure, setOfDOFS, groupElements, groupNodes);

    return surface;
}

void SetDBCPatchTestIGAIGA(NuTo::Structure &myStructure, int groupNodesSlave, int groupNodesMaster, int &countDBC)
{
    Eigen::Vector2d direction(0,0);
    double xMin(0.), xMax(0.);
    double yMin(0.), yMax(0.);
    auto LambdaNodes = [&](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double x = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[0];
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            bool xR = false;
            bool yR = false;

            if (x >= xMin - Tol && x <= xMax + Tol)
            {
                xR = true;
            }
            if (y >= yMin - Tol && y <= yMax + Tol)
            {
                yR = true;
            }

            if (xR == true && yR == true) return true;
        }
        return false;
    };

    // ===> master bottom <=== //
    xMin = 0.; xMax = std::numeric_limits<double>::infinity();
    yMin = yMax = -1.;
    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, LambdaNodes);
    direction << 0.,1.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesMasterLower, direction, 0.0);

    NuTo::FullVector<int,Eigen::Dynamic> rMembersMasterBottom;
    myStructure.NodeGroupGetMembers(groupNodesMasterLower, rMembersMasterBottom);
    direction << 1.,0.;
    countDBC = myStructure.ConstraintLinearSetDisplacementNode(rMembersMasterBottom(0), direction, 0.0);

    NuTo::FullVector<int, Eigen::Dynamic> members;

    // ===> PATCH TEST BOUNDARY <=== //
    double Stress = 10.;
    myStructure.SetNumLoadCases(1);
    // ===> slave top <=== //
    xMin = 1.; xMax = 2.;
    yMin = 1.; yMax = 1.;
    int groupNodesSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveUpper, LambdaNodes);
    members = myStructure.GroupGetMemberIds(groupNodesSlaveUpper);
    int groupElementsSlaveUpper = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveUpper, groupNodesSlaveUpper, false);
    members = myStructure.GroupGetMemberIds(groupElementsSlaveUpper);
    myStructure.LoadSurfacePressureCreate2D(0, groupElementsSlaveUpper, groupNodesSlaveUpper, Stress);

    countDBC++;
}

void MeshTyingIGA_IGA_Patch(const std::string &resultDir,
                            int rDegree,
                            double rPenalty,
                            NuTo::eIntegrationType rIntegrationType,
                            int rContactAlgo,
                            int refinementXSlave, int refinementYSlave,
                            int refinementElXMaster, int refinementElYMaster)
{
    NuTo::Structure myStructure(2);
    myStructure.SetNumTimeDerivatives(0);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure.SetNumProcessors(numThreads);
#endif

    /////////////////////////////////////////////////////////////////////
    // ====> create SLAVE mesh from gmsh (smth. impacting a rectangle) //
    /////////////////////////////////////////////////////////////////////
    int groupNodesSlave = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsSlave = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    double startX = 1.; double startY = 0.;
    double length = 1.; double height = 1.;
    buildRect2DIGA(rDegree, startX, startY, height, length,
                   groupElementsSlave, groupNodesSlave,  refinementXSlave, refinementYSlave,
                   &myStructure);

    // ===> build contact elements slave elements
    auto LambdaGetSlaveNodesLower = [](NuTo::NodeBase* rNodePtr) -> bool
    {
        double Tol = 1.e-6;
        if (rNodePtr->GetNum(NuTo::Node::eDof::COORDINATES)>0)
        {
            double y = rNodePtr->Get(NuTo::Node::eDof::COORDINATES)[1];
            if (y >= 0. - Tol && y <= 0. + Tol)
            {
                return true;
            }
        }
        return false;
    };  // LambdaGetSlaveNodesLower

    int groupNodesSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesSlaveLower, groupNodesSlave, LambdaGetSlaveNodesLower);

    NuTo::FullMatrix<int, Eigen::Dynamic> members = myStructure.GroupGetMemberIds(groupNodesSlaveLower);
    std::cout << "Members Nodes Slave lower: \n" << members << std::endl;

    int groupElementsSlaveLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsSlaveLower, groupNodesSlaveLower, false);
    members = myStructure.GroupGetMemberIds(groupElementsSlaveLower);
    std::cout << "Members Elements Slave lower: \n" << members << std::endl;

    //////////////////////////////
    // ====> create MASTER mesh //
    //////////////////////////////
    int groupNodesMaster = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    int groupElementsMaster = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    startX =  1.; startY = -1.;
    length =  1.; height =  1.;
    buildRect2DIGA(rDegree, startX, startY, height, length,
                   groupElementsMaster, groupNodesMaster,  refinementElXMaster,
                   refinementElYMaster, &myStructure);

    int groupNodesMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeFunction(groupNodesMasterLower, groupNodesMaster, LambdaGetSlaveNodesLower);

    members = myStructure.GroupGetMemberIds(groupNodesMasterLower);
    std::cout << "Members Nodes Master lower: \n" << members << std::endl;

    int groupElementsMasterLower = myStructure.GroupCreate(NuTo::eGroupId::Elements);
    myStructure.GroupAddElementsFromNodes(groupElementsMasterLower, groupNodesMasterLower, false);
    members = myStructure.GroupGetMemberIds(groupElementsMasterLower);
    std::cout << "Members Elements Master lower: \n" << members << std::endl;

    int countDBC;
    Eigen::VectorXd min1Coords(2,1);
    min1Coords << 0,0;
    Eigen::VectorXd max2Coords(2,1);
    max2Coords << 4,0;
    SetDBCPatchTestIGAIGA(myStructure, groupNodesSlave, groupNodesMaster, countDBC);

    ///////////////////
    // ===> material //
    ///////////////////

    double Thickness = 1.;
    int section = myStructure.SectionCreate("PLANE_STRESS");
    myStructure.SectionSetThickness(section, Thickness);

    double E   = 1.e5;
    double nue = 0.3;

    int constitutiveLaw = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    myStructure.ConstitutiveLawSetParameterDouble(constitutiveLaw, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, nue);
    myStructure.ElementTotalSetSection(section);
    myStructure.ElementTotalSetConstitutiveLaw(constitutiveLaw);

    ////////////////////////////
    // ===> CONTACT ELEMENTS  //
    ////////////////////////////

    int constitutiveLawPC = myStructure.ConstitutiveLawCreate("Contact_Constitutive_Law");
    std::function<double(double)> constitutiveContactLaw  =
    [rPenalty](double rGap) -> double
    {
        if(rGap<0)
            return rPenalty*rGap;
        else
            return 0.;
    };

    std::function<double(double)> constitutiveContactLawDerivative =
    [rPenalty](double rGap) -> double
    {
        if(rGap<=0)
            return rPenalty;
        else
            return 0.;
    };

    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC,
                                                    NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_FUNCTION,
                                                    constitutiveContactLaw);
    myStructure.ConstitutiveLawSetParameterFunction(constitutiveLawPC,
                                                    NuTo::Constitutive::eConstitutiveParameter::CONSTITUTIVE_LAW_DERIVATIVE_FUNCTION,
                                                    constitutiveContactLawDerivative);

    // *** mesh tying *** //

    int cceID = myStructure.ContactElementsCreate<2,2>(groupElementsSlaveLower, groupNodesSlaveLower,
                                                       groupElementsMasterLower, groupNodesMasterLower,
                                                       rIntegrationType, rContactAlgo, constitutiveLawPC);

    NuTo::ElementBase* elementBase = myStructure.ElementGetElementPtr(cceID);
    NuTo::ContinuumContactElement<2,2> &contactElement = elementBase->AsContinuumContactElement22();

    Eigen::MatrixXd D, M;
    std::unordered_map<int, int> mappingGlobal2LocalSlaveNode, mappingGlobal2LocalMasterNode;
    contactElement.ComputeMeshTyingMatrix(D, M, mappingGlobal2LocalSlaveNode, mappingGlobal2LocalMasterNode);

    std::cout << D << std::endl;

    std::cout << "===============================" << std::endl;

    std::cout << M << std::endl;

    std::vector<int> indicesBig;
    for(int j = 0; j < M.cols(); j++)
    {
        if(M.col(j).norm() > 1.e-10) indicesBig.push_back(j);
    }

    if(0)
    {
        for(int dof = 0; dof < myStructure.GetDimension(); dof++)
        {
            for(auto &itSlaveRow : mappingGlobal2LocalSlaveNode) // global
            {
                int localSlaveNodeIdRow  = itSlaveRow.second;

                std::unordered_map<int, int>::iterator itSlaveCol = mappingGlobal2LocalSlaveNode.begin();
                int localSlaveNodeIdCol = itSlaveCol->second;
                int globalSlaveNodeId = itSlaveCol->first;
                // define equation
                myStructure.ConstraintLinearEquationCreate(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, D(localSlaveNodeIdRow, localSlaveNodeIdCol), 0.);
                std::cout << D(localSlaveNodeIdRow, localSlaveNodeIdCol) << " ";
                itSlaveCol++;

                for(;itSlaveCol != mappingGlobal2LocalSlaveNode.end() ;itSlaveCol++)
                {
                    localSlaveNodeIdCol  = itSlaveCol->second;
                    globalSlaveNodeId = itSlaveCol->first;
                    // add slave ingredients
                    myStructure.ConstraintLinearEquationAddTerm(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, D(localSlaveNodeIdRow, localSlaveNodeIdCol));
                    std::cout << D(localSlaveNodeIdRow, localSlaveNodeIdCol) << " ";
                    // master side
                }

                std::cout << " || ";
                for(auto &itMaster : mappingGlobal2LocalMasterNode)
                {
                    int globalMasterNodeId = itMaster.first;
                    int localMasterNodeIdCol  = itMaster.second;
                    // add master ingredients
                    for(auto &it: indicesBig)
                        if(it == localMasterNodeIdCol)
                            myStructure.ConstraintLinearEquationAddTerm(countDBC, globalMasterNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, M(localSlaveNodeIdRow, localMasterNodeIdCol));

                    std::cout << M(localSlaveNodeIdRow, localMasterNodeIdCol) << " ";
                }
                std::cout << std::endl;

                countDBC++;
            }
        }
    }
    else
    {
        Eigen::MatrixXd D_small, M_small;
        D_small.resize(mappingGlobal2LocalSlaveNode.size(), mappingGlobal2LocalSlaveNode.size());
        M_small.resize(mappingGlobal2LocalSlaveNode.size(), mappingGlobal2LocalMasterNode.size());
        int i = 0, jSlave = 0, jMaster  = 0;
        for(auto &itSlaveRow : mappingGlobal2LocalSlaveNode)
        {
            int row = itSlaveRow.second;
            jSlave = 0;
            for(auto &itSlaveCol : mappingGlobal2LocalSlaveNode)
            {
                int col = itSlaveCol.second;
                D_small(i,jSlave) = D(row,col);
                jSlave++;
            }
            jMaster = 0;
            for(auto &itMasterCol : mappingGlobal2LocalMasterNode)
            {
                int colMaster = itMasterCol.second;
                M_small(i,jMaster) = M(row,colMaster);
                jMaster++;
            }
            i++;
        }

        std::vector<int> indices;
        for(int j = 0; j < M_small.cols(); j++)
        {
            if(M_small.col(j).norm() > 1.e-10) indices.push_back(j);
        }

        std::cout << "D_small:\n";
        std::cout << std::endl << D_small << std::endl << std::endl;
        std::cout << "M_small:\n";
        std::cout << M_small << std::endl << std::endl;
        std::cout << "Ml:\n";
        std::cout << M << std::endl << std::endl;
        Eigen::MatrixXd P = (D_small.inverse())*M_small;
        std::cout << D_small.inverse()*(D_small) << std::endl << std::endl;
        std::cout << "P:\n";
        std::cout << P << std::endl << std::endl;

        i = 0;
        for(auto &itSlaveRow : mappingGlobal2LocalSlaveNode)
        {
            int globalSlaveNodeId = itSlaveRow.first;
            for(int dof = 0; dof < myStructure.GetDimension(); dof++)
            {
                myStructure.ConstraintLinearEquationCreate(countDBC, globalSlaveNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, 1., 0.);
                int j = 0;
                for(auto &itMasterCol : mappingGlobal2LocalMasterNode)
                {
                    for(auto &it: indices)
                    {
                        if(it == j)
                        {
                            std::cout << j << std::endl;
                            int globalMasterNodeId = itMasterCol.first;
                            myStructure.ConstraintLinearEquationAddTerm(countDBC, globalMasterNodeId, NuTo::Node::eDof::DISPLACEMENTS, dof, P(i, j));
                        }
                    }
                    j++;
                }
                countDBC++;
            }
            i++;
        }
    }

    myStructure.ElementDelete(cceID);

    ///////////////////
    // ===> Solution //
    ///////////////////
    //set result directory
    if (boost::filesystem::exists(resultDir))
    {
        if (boost::filesystem::is_directory(resultDir))
        {
            boost::filesystem::remove_all(resultDir);
        }
    }

    // create result directory
    boost::filesystem::create_directory(resultDir);

    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    double timeStep = 1.;
    double simulationTime = 1.;

    myIntegrationScheme.SetResultDirectory(resultDir, true);

    myIntegrationScheme.SetMinTimeStepPlot(1.);
    myIntegrationScheme.SetLastTimePlot(0.);

    myIntegrationScheme.SetToleranceForce(1.e-10);
    myIntegrationScheme.SetMaxNumIterations(50);
    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.Solve(simulationTime);

#ifdef ENABLE_VISUALIZE
    int visualizationGroup = myStructure.GroupUnion(groupElementsSlave, groupElementsMaster);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    myStructure.ElementGroupExportVtkDataFile(visualizationGroup, resultDir+"/Elements.vtu", true);
#endif
}

int main()
{
    std::string resultDir = "./ResultsStatic_Tying_IGA_IGA_Patch_Mortar_12IP";
    int contactAlgorithm = 0;
    int degree = 2;
    double penalty = 1.e9;

    MeshTyingIGA_IGA_Patch(resultDir,
                           degree,
                           penalty,
                           NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 1, 1, 3, 1);

//    resultDir = "./ResultsStatic_Tying_IGA_IGA_Patch_NonMortar";
//    contactAlgorithm = 1;
//    MeshTyingIGA_IGA_Patch(resultDir,
//                             degree,
//                             penalty,
//                             NuTo::eIntegrationType::IntegrationType1D2NGauss12Ip, contactAlgorithm, 1, 1, 2, 1);

    return 0;
}
