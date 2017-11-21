#define TRILINOS_USE_INSTALLED_LIBS
//#undef TRILINOS_USE_INSTALLED_LIBS

#include <mpi.h>

#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <map>
//#include <eigen3/Eigen/Core>
//#include <eigen3/Eigen/Dense>

#include "json.hpp"
#include "base/Exception.h"

#ifdef TRILINOS_USE_INSTALLED_LIBS

#include <Epetra_MpiComm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Export.h>
#include <Epetra_LinearProblem.h>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Export.hpp>
#include <Amesos2.hpp>
#include <Amesos2_Meta.hpp>
//#include <Amesos2_TpetraCrsMatrix_MatrixAdapter.hpp>
//#include <Amesos2_TpetraMultiVecAdapter.hpp>

#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>


#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>
#include <Amesos_Mumps.h>
#include <Amesos_ConfigDefs.h>
#include <Teuchos_ParameterList.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack.h>
#include <Ifpack_AdditiveSchwarz.h>

#else
#include <trilinos/Epetra_MpiComm.h>
#include <trilinos/Epetra_CrsMatrix.h>
#include <trilinos/Epetra_Map.h>
#include <trilinos/Epetra_Export.h>
#include <trilinos/Epetra_LinearProblem.h>

#include <trilinos/Tpetra_DefaultPlatform.hpp>
#include <trilinos/Tpetra_CrsMatrix.hpp>
#include <trilinos/Tpetra_CrsGraph.hpp>
#include <trilinos/Tpetra_Map.hpp>
#include <trilinos/Tpetra_Export.hpp>
#include <trilinos/Amesos2.hpp>
#include <trilinos/BelosTpetraAdapter.hpp>
#include <trilinos/BelosSolverFactory.hpp>

#include <trilinos/Teuchos_RCP.hpp>
#include <trilinos/Teuchos_Comm.hpp>

#include <trilinos/AztecOO.h>
#include <trilinos/Amesos.h>
#include <trilinos/Amesos_BaseSolver.h>
#include <trilinos/Amesos_Mumps.h>
#include <trilinos/Amesos_ConfigDefs.h>
#include <trilinos/Teuchos_ParameterList.hpp>
#include <trilinos/BelosLinearProblem.hpp>
#include <trilinos/BelosBlockGmresSolMgr.hpp>
#include <trilinos/BelosPseudoBlockGmresSolMgr.hpp>
#include <trilinos/BelosEpetraAdapter.hpp>
#include <trilinos/Ifpack.h>
#include <trilinos/Ifpack_AdditiveSchwarz.h>
#endif

#include "mechanics/groups/Group.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/feti/StructureFeti.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/AdditiveOutput.h"

#include "../applications/custom/TestClasses/ConversionTools.h"
#include "../applications/custom/TestClasses/StructureMesh.h"
#include "../applications/custom/TestClasses/MeshFileGenerator.h"


using NuTo::Interpolation::eShapeType;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Node::eDof;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::eDirection;

using Teuchos::RCP;
using Teuchos::rcp;

void run_mesh_test(Epetra_MpiComm rComm, std::string rFileName)
{
    //******************************************
    //*            get MPI setting             *
    //******************************************
    int numProc = rComm.NumProc();
    int rank = rComm.MyPID();

    RCP<const Teuchos::Comm<int>> commTeuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    numProc = commTeuchos->getSize();
    rank = commTeuchos->getRank();

    //******************************************
    //*       set structure constants          *
    //******************************************
    int dim = 2;
    const double thickness = 1.0;
    const double lengthX = 60.;
    const double lengthY = 10.;
    const double youngsModulus = 2.1e5;
    const double poissonsRatio = 0.3;
    double displacementValue = 0.0;
    double forceValue = 1.e6;
    NuTo::Node::eDof dofType = NuTo::Node::eDof::DISPLACEMENTS;


    //******************************************
    //*           set up structure             *
    //******************************************
    StructureMesh structure(dim);
//    std::string meshFile = "../meshes/jsonMesh_Test_NodeIDs" + std::to_string(rank) + ".mesh";
    std::string meshFile = rFileName + "_" + std::to_string(rank+1);


    //******************************************
    //*       set interpolation type           *
    //******************************************
//    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
//    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
//    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);


    //******************************************
    //*           import geometry              *
    //******************************************
//    structure.importMyMeshJson(meshFile, interpolationTypeId);
    structure.importMyMeshJson(meshFile);


    //******************************************
    //*         set constitutive law           *
    //******************************************
    const int linElasticStress = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(linElasticStress, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(linElasticStress, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);

    const int heatConduct = structure.ConstitutiveLawCreate(eConstitutiveType::HEAT_CONDUCTION);
    structure.ConstitutiveLawSetParameterDouble(heatConduct, eConstitutiveParameter::HEAT_CAPACITY, 1e-3);

    const int thermStrains = structure.ConstitutiveLawCreate(eConstitutiveType::THERMAL_STRAINS);
    structure.ConstitutiveLawSetParameterDouble(thermStrains, NuTo::Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 15e-6);

    structure.ElementTotalSetConstitutiveLaw(linElasticStress);

//    int additive_input_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_INPUT_EXPLICIT);
//    int additive_output_id = structure.ConstitutiveLawCreate(eConstitutiveType::ADDITIVE_OUTPUT);

//    NuTo::AdditiveInputExplicit* additive_input = static_cast<NuTo::AdditiveInputExplicit*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_input_id));
//    NuTo::AdditiveOutput* additive_output = static_cast<NuTo::AdditiveOutput*>(structure.ConstitutiveLawGetConstitutiveLawPtr(additive_output_id));
//    NuTo::ConstitutiveBase* lin_elastic = structure.ConstitutiveLawGetConstitutiveLawPtr(linElasticStress);
//    NuTo::ConstitutiveBase* thermal_strains = structure.ConstitutiveLawGetConstitutiveLawPtr(thermStrains);
//    NuTo::ConstitutiveBase* heat_conduction = structure.ConstitutiveLawGetConstitutiveLawPtr(heatConduct);

//    additive_input->AddConstitutiveLaw(*lin_elastic);
//    additive_input->AddConstitutiveLaw(*thermal_strains, NuTo::Constitutive::eInput::ENGINEERING_STRAIN);

//    additive_output->AddConstitutiveLaw(*additive_input);
//    additive_output->AddConstitutiveLaw(*heat_conduction);

//    structure.ElementTotalSetConstitutiveLaw(additive_output_id);
//    structure.ElementGroupSetConstitutiveLaw(group, additive_output_id);


    //******************************************
    //*            define section              *
    //******************************************
    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);


    //******************************************
    //*         define constraints             *
    //******************************************
    auto& groupNodesLeftBoundary = structure.GroupGetNodesAtCoordinate(eDirection::X, 0.);
    structure.Constraints().Add(eDof::DISPLACEMENTS, NuTo::Constraint::Component(groupNodesLeftBoundary, {eDirection::X, eDirection::Y}, displacementValue));
//    structure.Constraints().Add(eDof::DISPLACEMENTS, NuTo::Constraint::Component(nodesMiddleLowerBoundary, {eDirection::Y}, displValue));
//    structure.Constraints().Add(eDof::DISPLACEMENTS, NuTo::Constraint::Component(groupNodesLeftBoundary, eDirection::X, 0.0));
//    structure.Constraints().Add(eDof::TEMPERATURE, NuTo::Constraint::Value(groupNodesLeftBoundary));


    //******************************************
    //*            define loads                *
    //******************************************
//    auto& groupNodesRightBoundary = structure.GroupGetNodesAtCoordinate(eDirection::X, lengthX);
//    auto& nodesMiddleLowerBoundary = structure.GroupGetNodeRadiusRange(Eigen::Vector2d(lengthX/2., 0));
//    auto& nodesMiddleUpperBoundary = structure.GroupGetNodeRadiusRange(Eigen::Vector2d(lengthX/2., lengthY));
    auto& nodesRightUpperBoundary = structure.GroupGetNodeRadiusRange(Eigen::Vector2d(lengthX, lengthY));
    auto& nodesRightLowerBoundary = structure.GroupGetNodeRadiusRange(Eigen::Vector2d(lengthX, 0));
//    structure.LoadCreateNodeGroupForce(&groupNodesRightBoundary, Eigen::Vector2d::UnitX(), forceValue);
//    structure.LoadCreateNodeGroupForce(&groupNodesRightBoundary, Eigen::Vector2d(1, 0.1), forceValue);
//    structure.LoadCreateNodeGroupForce(&nodesMiddleLowerBoundary, -Eigen::Vector2d::UnitY(), 0.5*forceValue);
//    structure.LoadCreateNodeGroupForce(&nodesMiddleLowerBoundary, -Eigen::Vector2d::UnitY(), 2.5*forceValue);
//    structure.LoadCreateNodeGroupForce(&nodesMiddleUpperBoundary, Eigen::Vector2d::UnitY(), 2.5*forceValue);
    structure.LoadCreateNodeGroupForce(&nodesRightUpperBoundary, Eigen::Vector2d(0.2,1), forceValue);
    structure.LoadCreateNodeGroupForce(&nodesRightLowerBoundary, Eigen::Vector2d(0.2,-1), forceValue);


    //******************************************
    //*        compute active system           *
    //******************************************
    Eigen::SparseMatrix<double> A_JJ = structure.BuildGlobalHessian0().JJ.ExportToEigenSparseMatrix();
    Eigen::MatrixXd r_J = (structure.BuildGlobalInternalGradient() - structure.BuildGlobalExternalLoadVector()).J.Export();
    structure.SetVerboseLevel(10);
    structure.Info();


    //******************************************
    //*        generate dof setting            *
    //******************************************
    structure.generateNodeToDofMapping();
    structure.generateDofClassification();
//    structure.gatherNodeToDofMapping_allProcesses(numProc);
    structure.gatherNodeToDofMapping(numProc, rank);


    //******************************************
    //*          get dof mappings              *
    //******************************************
//    std::vector<std::map<NuTo::Node::eDof, std::map<int, int>>> local2GlobalDof = structure.getLocalToGlobalDofMapping();
    std::map<NuTo::Node::eDof, std::map<int, int>> myLocal2GlobalDofs = structure.getMyLocalToGlobalDofMapping();
    std::map<NuTo::Node::eDof, std::vector<int>> myMasterGlobalDofs = structure.getMyMasterGlobalDofIDs();
    std::map<NuTo::Node::eDof, std::vector<int>> myMasterGlobalActiveDofs = structure.getMyMasterGlobalActiveDofIDs();


    //******************************************
    //*     create overlapping index map       *
    //******************************************
//    std::vector<int> myGlobalActiveDofIDs = structure.map2Vector(local2GlobalActiveDof[rank][dofType]);
    std::vector<int> myGlobalDofIDs = structure.map2Vector(myLocal2GlobalDofs[dofType]);
    int* myGlobalDofIDs_arr = &myGlobalDofIDs[0];
    Epetra_Map overlappingMap(-1, myGlobalDofIDs.size(), myGlobalDofIDs_arr, 0, rComm);
    RCP<Tpetra::Map<int, int>> overlappingMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myGlobalDofIDs_arr, myGlobalDofIDs.size(), 0, commTeuchos));
//    std::cout << overlappingMap_tpetra.description() << std::endl;


    //******************************************
    //*       create owning index map          *
    //******************************************
//    std::vector<int> myOwningGlobalActiveDofIDs = masterGlobalActiveDofs[rank][dofType];
    std::vector<int> myOwningGlobalActiveDofIDs = myMasterGlobalActiveDofs[dofType];
    int* myOwningGlobalActiveDofIDs_arr = &myOwningGlobalActiveDofIDs[0];
    Epetra_Map owningMap(-1, myOwningGlobalActiveDofIDs.size(), myOwningGlobalActiveDofIDs_arr, 0, rComm);
    RCP<Tpetra::Map<int, int>> owningMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myOwningGlobalActiveDofIDs_arr, myOwningGlobalActiveDofIDs.size(), 0, commTeuchos));


    owningMap.Print(std::cout);
    overlappingMap.Print(std::cout);

    //******************************************
    //*         create index graphs            *
    //******************************************
    int maxNonZeros = 18;   //got by interpolation type (order), e.g. QUAD2 + EQUIDISTANT1 => 18
    Epetra_CrsGraph owningGraph(Epetra_DataAccess::Copy, owningMap, maxNonZeros, false);
    Epetra_CrsGraph overlappingGraph(Epetra_DataAccess::Copy, overlappingMap, maxNonZeros, false);
    RCP<Tpetra::CrsGraph<int, int>> owningGraph_tpetra = rcp(new Tpetra::CrsGraph<int,int>(owningMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
    RCP<Tpetra::CrsGraph<int, int>> overlappingGraph_tpetra = rcp(new Tpetra::CrsGraph<int, int>(overlappingMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
    int errCode = 0;    //TODO: error handling
    for (int k=0; k<A_JJ.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A_JJ,k); it; ++it)
        {
            // describe position of entries
            errCode = overlappingGraph.InsertGlobalIndices(myGlobalDofIDs_arr[it.row()], 1, &myGlobalDofIDs_arr[it.col()]);
            overlappingGraph_tpetra->insertGlobalIndices(myGlobalDofIDs_arr[it.row()], 1, &myGlobalDofIDs_arr[it.col()]);
//            int col = it.col();
//            errCode = overlappingGraph.InsertMyIndices(it.row(), 1, &col);
        }
    }

    //******************************************
    //*   define inter-process communication   *
    //*      for local-to-global indices       *
    //******************************************
    Epetra_Export exporter(overlappingMap, owningMap);
    errCode = overlappingGraph.FillComplete();
    errCode = owningGraph.Export(overlappingGraph, exporter, Insert); //inter-process communication of matrix structure
    errCode = owningGraph.FillComplete();

    RCP<Tpetra::Export<int, int>> exporter_tpetra = rcp(new Tpetra::Export<int, int>(overlappingMap_tpetra, owningMap_tpetra));
    overlappingGraph_tpetra->fillComplete();
    owningGraph_tpetra->doExport(*overlappingGraph_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    owningGraph_tpetra->fillComplete();


    //******************************************
    //*  initialize Trilinos matrix and vector *
    //******************************************
    Epetra_CrsMatrix globalA_JJ(Epetra_DataAccess::Copy, owningGraph);
    Epetra_Vector globalRhsVector(owningMap);
    globalRhsVector.PutScalar(0.0);
    Tpetra::CrsMatrix<double, int, int> globalA_JJ_tpetra(owningGraph_tpetra);
    Tpetra::Vector<double, int, int> globalRhsVector_tpetra(owningMap_tpetra);
    globalRhsVector_tpetra.putScalar(0.0);


    //******************************************
    //*    conversion from NuTo to Trilinos    *
    //******************************************
    ConversionTools converter(rComm);
    Eigen::SparseMatrix<double, Eigen::RowMajor> A_JJ_rowMajor(A_JJ);
    Epetra_CrsMatrix localA_JJ = converter.convertEigen2EpetraCrsMatrix(A_JJ_rowMajor, overlappingGraph, true);
    Epetra_Vector localRhsVector = converter.convertEigen2EpetraVector(r_J, overlappingMap);
    errCode = globalA_JJ.Export(localA_JJ, exporter, Add);
    errCode = globalRhsVector.Export(localRhsVector, exporter, Insert);
    globalRhsVector.Scale(-1.);

    ConversionTools converter2;
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> localA_JJ_tpetra = converter2.convertEigen2TpetraCrsMatrix(A_JJ_rowMajor, overlappingGraph_tpetra, true);
    Teuchos::RCP<Tpetra::Vector<double, int, int>> localRhsVector_tpetra = converter2.convertEigen2TpetraVector(r_J, overlappingMap_tpetra);
    globalA_JJ_tpetra.doExport(*localA_JJ_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::ADD);
    globalRhsVector_tpetra.doExport(*localRhsVector_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    globalRhsVector_tpetra.scale(-1.);


    //******************************************
    //*        solve complete problem          *
    //******************************************
    Epetra_MultiVector sol = structure.solveSystem(globalA_JJ, globalRhsVector, false);
    sol.Print(std::cout);   //print solution


    //******************************************
    //*          visualize solution            *
    //******************************************
    std::vector<double> sol_std = converter.convertEpetraMultiVector2StdVector(sol, 0, false);
    std::vector<int> myOwningGlobalDofIDs = myMasterGlobalDofs[dofType];
//    myOwningGlobalDofIDs = myGlobalDofIDs;
    std::vector<std::vector<int>> node2Dofs_vector = structure.map2Vector(structure.getMyNodeDOFs(dofType));
    std::map<NuTo::Node::eDof, std::vector<int>> masterGlobalDependentDofs = structure.getMyMasterGlobalDependentDofIDs();
    // add previously defined dependent dof values
    for (int j = 0; j < masterGlobalDependentDofs[dofType].size(); ++j)
    {
        sol_std.push_back(displacementValue);
    }

    for (int mogi : myOwningGlobalDofIDs)
        std::cout << rank << ": " << mogi << std::endl;

    structure.visualizeSerializedParticularSolution(sol_std, myOwningGlobalDofIDs, node2Dofs_vector, "result" + std::to_string(rank) + ".vtu", numProc);
//    structure.visualizeSolution(sol_std, "result_default" + std::to_string(rank) + ".vtu");
}


Teuchos::RCP<Tpetra::MultiVector<double, int, int>> solveSystem_tpetra(Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int>> rA, Teuchos::RCP<const Tpetra::Vector<double, int, int>> rRhs, Teuchos::RCP<Tpetra::Vector<double, int, int>> rLhs, bool rIterative = false)
{
    if (rIterative)
    {
//        using Teuchos::ParameterList;

//        // Make an empty new parameter list.
//        RCP<ParameterList> solverParams = Teuchos::parameterList();

//        // Set some GMRES parameters.
//        solverParams->set ("Num Blocks", 40);       //Maximum number of Krylov vectors to store, also restart length
//        solverParams->set ("Maximum Iterations", 400);
//        solverParams->set ("Convergence Tolerance", 1.0e-8);

//        // Create the GMRES solver using a "factory" and
//        // the list of solver parameters created above.
//        Belos::SolverFactory<Tpetra::Vector<>::scalar_type, Tpetra::Vector<double, int, int>, Tpetra::CrsMatrix<double, int, int>> factory;
//        RCP<Belos::SolverManager<Tpetra::Vector<>::scalar_type, Tpetra::Vector<double, int, int>, Tpetra::CrsMatrix<double, int, int>>> solver = factory.create ("GMRES", solverParams);

//        // Create a LinearProblem struct with the problem to solve.
//        // A, X, B, and M are passed by (smart) pointer, not copied.
//        typedef Belos::LinearProblem<Tpetra::Vector<>::scalar_type, Tpetra::Vector<double, int, int>, Tpetra::CrsMatrix<double, int, int>> BelosLinearProblem;
//        RCP<BelosLinearProblem> problem = rcp (new BelosLinearProblem (rA, rLhs, rRhs));

//        // You don't have to call this if you don't have a preconditioner.
//        // If M is null, then Belos won't use a (right) preconditioner.
//        //problem->setRightPrec (M);

//        // Tell the LinearProblem to make itself ready to solve.
//        problem->setProblem ();

//        // Tell the solver what problem you want to solve.
//        solver->setProblem (problem);

//        Belos::ReturnType result = solver->solve();

//        if (result == Belos::Converged)
//        {
//            Tpetra::Vector<double, int, int> result =  *solver->getProblem().getLHS().get();
//            return rLhs;
//        }
//        else
//        {
//            return rcp(new Tpetra::MultiVector<double, int, int>());
//        }

    }
    else
    {
        // Before we do anything, check that SuperLU is enabled
        if( !Amesos2::query("SuperLU") ){
          std::cerr << "SuperLU not enabled.  Exiting..." << std::endl;
          return rcp(new Tpetra::MultiVector<double, int, int>());      // Otherwise CTest will pick it up as
                                                                        // failure, which it isn't really
        }

        typedef Tpetra::Map<>::local_ordinal_type loc_ord_type;
        typedef Tpetra::Map<>::global_ordinal_type glob_ord_type;
        typedef Tpetra::CrsMatrix<double, loc_ord_type, glob_ord_type> mat;
        typedef Tpetra::MultiVector<double, loc_ord_type, glob_ord_type> vec;

        // Create solver interface to Superlu with Amesos2 factory method
        RCP<Amesos2::Solver<mat, vec>> solver = Amesos2::create<mat, vec>("SuperLU", rA, rLhs, rRhs);
        solver->symbolicFactorization().numericFactorization().solve();

        std::ostream &out = std::cout;
        RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
        *fos << "Solution :" << std::endl;
        rLhs->describe(*fos,Teuchos::VERB_EXTREME);
        *fos << std::endl;

        return rLhs;
    }

    return rLhs;
}

Teuchos::RCP<Tpetra::MultiVector<double, int, int>> solveSystem_tpetra(Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int>> rA, Teuchos::RCP<const Tpetra::Vector<double, int, int>> rRhs, bool rIterative = false)
{
    Teuchos::RCP<Tpetra::Vector<double, int, int>> lhs = rcp(new Tpetra::Vector<double, int, int>(*rRhs.get()));

    return solveSystem_tpetra(rA, rRhs, lhs, rIterative);
}


void run_generation_test()
{
    std::string gmshFileName = "../meshes/origMesh_20.msh_000001";
    gmshFileName = "../meshes/origMesh_200.msh";
    std::string meshFileName = "../meshes/origMesh_200";

    MeshFileGenerator gen;
//    gen.generateMeshFromGmsh(gmshFileName, meshFileName);
    std::vector<MeshFileGenerator::NuTo_DofTypes> dofTypes;
    MeshFileGenerator::NuTo_InterpolationOrders interpolationOrder = MeshFileGenerator::NuTo_InterpolationOrders::EQUIDISTANT1;
    gen.generateMeshFilesFromGmsh(gmshFileName, meshFileName, dofTypes, interpolationOrder);
}


int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    Epetra_MpiComm comm(MPI_COMM_WORLD);

//    run_generation_test();

    std::string fileNames = "../meshes/origMesh_200.mff";
//    fileNames = "../meshes/jsonMesh_Test_NodeIDs.mesh";
    run_mesh_test(comm, fileNames);



    MPI_Finalize();
    return 0;
}































