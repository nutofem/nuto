#define TRILINOS_USE_INSTALLED_LIBS
//#undef TRILINOS_USE_INSTALLED_LIBS

#include <mpi.h>

#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include <map>

#include "json.hpp"

#ifdef TRILINOS_USE_INSTALLED_LIBS

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixFactory.hpp>

#include <Amesos2.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <Ifpack.h>
#include <Ifpack_AdditiveSchwarz.h>

#include <zoltan.h>

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


#include "base/Group.h"

#include "mechanics/mesh/MeshFem.h"
#include "mechanics/mesh/MeshFemDofConvert.h"

#include "mechanics/interpolation/InterpolationQuadLinear.h"
#include "mechanics/interpolation/InterpolationTrussLinear.h"
#include "mechanics/integrationtypes/IntegrationTypeTensorProduct.h"

#include "mechanics/constitutive/laws/LinearElastic.h"

#include "mechanics/integrands/MomentumBalance.h"
#include "mechanics/integrands/NeumannBc.h"

#include "mechanics/cell/Cell.h"
#include "mechanics/cell/SimpleAssember.h"

#include "../TestClasses/ConversionTools.h"
#include "../TestClasses/StructureMesh.h"
#include "../TestClasses/MeshFileGenerator.h"
#include "../TestClasses/ZoltanMesh.h"


using NuTo::Interpolation::eShapeType;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Node::eDof;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::eDirection;

using Teuchos::RCP;
using Teuchos::rcp;


Teuchos::RCP<Tpetra::MultiVector<double, int, int>> solveSystem_tpetra(Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int>> rA, Teuchos::RCP<const Tpetra::Vector<double, int, int>> rRhs, Teuchos::RCP<Tpetra::Vector<double, int, int>> rLhs, bool rIterative = false, bool printSolution = true)
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
        // Before we do anything, check that solver is available
        //
        // Solver:  - KLU2
        //          - Mumps
        //          - SuperLU

        std::string solverType = "KLU2";
        if( !Amesos2::query(solverType) ){
          std::cerr << solverType << " not enabled.  Exiting..." << std::endl;
          return rcp(new Tpetra::MultiVector<double, int, int>());      // Otherwise CTest will pick it up as
                                                                        // failure, which it isn't really
        }

        typedef Tpetra::Map<>::local_ordinal_type loc_ord_type;
        typedef Tpetra::Map<>::global_ordinal_type glob_ord_type;
//        typedef Tpetra::CrsMatrix<double, loc_ord_type, glob_ord_type> mat;
//        typedef Tpetra::MultiVector<double, loc_ord_type, glob_ord_type> vec;
        typedef Tpetra::CrsMatrix<double, int, int> mat;
        typedef Tpetra::MultiVector<double, int, int> vec;

        // Create solver interface to Superlu with Amesos2 factory method
        RCP<Amesos2::Solver<mat, vec>> solver = Amesos2::create<mat, vec>(solverType, rA, rLhs, rRhs);
        solver->symbolicFactorization();
        solver->numericFactorization();
        solver->solve();
//        solver->symbolicFactorization().numericFactorization().solve();

        if (printSolution)
        {
            std::ostream &out = std::cout;
            RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
            *fos << "Solution :" << std::endl;
            rLhs->describe(*fos,Teuchos::VERB_EXTREME);
            *fos << std::endl;
        }

        return rLhs;
    }

    return rLhs;
}

Teuchos::RCP<Tpetra::MultiVector<double, int, int>> solveSystem_tpetra(Teuchos::RCP<const Tpetra::CrsMatrix<double, int, int>> rA, Teuchos::RCP<const Tpetra::Vector<double, int, int>> rRhs, bool rIterative = false)
{
//    Teuchos::RCP<Tpetra::Vector<double, int, int>> lhs = rcp(new Tpetra::Vector<double, int, int>(*rRhs.get()));
    Teuchos::RCP<Tpetra::Vector<double, int, int>> lhs = rcp(new Tpetra::Vector<double, int, int>(rA->getDomainMap()));

    return solveSystem_tpetra(rA, rRhs, lhs, rIterative);
}

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
//    RCP<Xpetra::Map<int, int>> overlappingMap_xpetra = rcp(new Xpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myGlobalDofIDs_arr, myGlobalDofIDs.size(), 0, commTeuchos));
//    std::cout << overlappingMap_tpetra.description() << std::endl;


    //******************************************
    //*       create owning index map          *
    //******************************************
//    std::vector<int> myOwningGlobalActiveDofIDs = masterGlobalActiveDofs[rank][dofType];
    std::vector<int> myOwningGlobalActiveDofIDs = myMasterGlobalActiveDofs[dofType];
    int* myOwningGlobalActiveDofIDs_arr = &myOwningGlobalActiveDofIDs[0];
    Epetra_Map owningMap(-1, myOwningGlobalActiveDofIDs.size(), myOwningGlobalActiveDofIDs_arr, 0, rComm);
    RCP<Tpetra::Map<int, int>> owningMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myOwningGlobalActiveDofIDs_arr, myOwningGlobalActiveDofIDs.size(), 0, commTeuchos));
//    RCP<Xpetra::Map<int, int>> owningMap_xpetra = rcp(new Xpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myOwningGlobalActiveDofIDs_arr, myOwningGlobalActiveDofIDs.size(), 0, commTeuchos));


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
//    RCP<Xpetra::CrsGraph<int, int>> owningGraph_xpetra = rcp(new Xpetra::CrsGraph<int,int>(owningMap_xpetra, maxNonZeros, Xpetra::ProfileType::DynamicProfile));
//    RCP<Xpetra::CrsGraph<int, int>> overlappingGraph_xpetra = rcp(new Xpetra::CrsGraph<int, int>(overlappingMap_xpetra, maxNonZeros, Xpetra::ProfileType::DynamicProfile));
    std::vector<int> columnIndices;
    int rowIndex = 0;
    int errCode = 0;    //TODO: error handling
    for (int k=0; k<A_JJ.outerSize(); ++k)
    {
        columnIndices.clear();
        for (Eigen::SparseMatrix<double>::InnerIterator it(A_JJ,k); it; ++it)
        {
            // describe position of entries
            errCode = overlappingGraph.InsertGlobalIndices(myGlobalDofIDs_arr[it.row()], 1, &myGlobalDofIDs_arr[it.col()]);
            overlappingGraph_tpetra->insertGlobalIndices(myGlobalDofIDs_arr[it.row()], 1, &myGlobalDofIDs_arr[it.col()]);
//            columnIndices.push_back(myGlobalDofIDs_arr[it.col()]);
//            rowIndex = myGlobalDofIDs_arr[it.row()];
//            int col = it.col();
//            errCode = overlappingGraph.InsertMyIndices(it.row(), 1, &col);
        }
//        Teuchos::ArrayView<int> colInds = Teuchos::arrayViewFromVector(columnIndices);
//        overlappingGraph_xpetra->insertGlobalIndices(rowIndex, colInds);
    }

    //******************************************
    //*   define inter-process communication   *
    //*      for local-to-global indices       *
    //******************************************
    Epetra_Export exporter(overlappingMap, owningMap);
    errCode = overlappingGraph.FillComplete();
    errCode = owningGraph.Export(overlappingGraph, exporter, Insert); //inter-process communication of matrix structure
    errCode = owningGraph.FillComplete();

    RCP<const Tpetra::Export<int, int>> exporter_tpetra = rcp(new Tpetra::Export<int, int>(overlappingMap_tpetra, owningMap_tpetra));
    overlappingGraph_tpetra->fillComplete();
    owningGraph_tpetra->doExport(*overlappingGraph_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    owningGraph_tpetra->fillComplete();

//    RCP<Xpetra::Export<int, int>> exporter_xpetra = rcp(new Xpetra::Export<int, int>(overlappingMap_xpetra, owningMap_xpetra));
//    overlappingGraph_xpetra->fillComplete();
//    owningGraph_xpetra->doExport(*overlappingGraph_xpetra.get(), *exporter_xpetra.get(), Xpetra::CombineMode::INSERT);
//    owningGraph_xpetra->fillComplete();


    //******************************************
    //*  initialize Trilinos matrix and vector *
    //******************************************
    Epetra_CrsMatrix globalA_JJ(Epetra_DataAccess::Copy, owningGraph);
    Epetra_Vector globalRhsVector(owningMap);
    globalRhsVector.PutScalar(0.0);
//    Tpetra::CrsMatrix<double, int, int> globalA_JJ_tpetra(owningGraph_tpetra);
//    Tpetra::Vector<double, int, int> globalRhsVector_tpetra(owningMap_tpetra);
//    globalRhsVector_tpetra.putScalar(0.0);
    RCP<Tpetra::CrsMatrix<double, int, int>> globalA_JJ_tpetra = rcp(new Tpetra::CrsMatrix<double, int, int>(owningGraph_tpetra));
    RCP<Tpetra::Vector<double, int, int>> globalRhsVector_tpetra = rcp(new Tpetra::Vector<double, int, int>(owningMap_tpetra));
    globalRhsVector_tpetra->putScalar(0.0);

//    RCP<Xpetra::CrsMatrix<double, int, int>> globalA_JJ_xpetra = Xpetra::CrsMatrixFactory<double, int, int, Kokkos::Compat::KokkosSerialWrapperNode>::Build(owningGraph_tpetra);
//    RCP<Xpetra::CrsMatrix<double, int, int>> globalA_JJ_xpetra = rcp(new Xpetra::CrsMatrix<double, int, int>(owningGraph_xpetra));
//    RCP<Xpetra::Vector<double, int, int>> globalRhsVector_xpetra = rcp(new Xpetra::Vectro<double, int, int>(owningMap_xpetra));
//    globalRhsVector_xpetra->putScalar(0.0);


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
//    RCP<Xpetra::CrsMatrix<double, int, int>> localA_JJ_xpetra = converter2.convertEigen2XpetraCrsMatrix<double, int, int>(A_JJ_rowMajor, overlappingGraph_tpetra, true);
//    RCP<Xpetra::Vector<double, int, int>> localRhsVector_xpetra = converter2.convertEigen2XpetraVector<double, int, int>(r_J, overlappingMap_tpetra);
    globalA_JJ_tpetra->doExport(*localA_JJ_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::ADD);
    globalRhsVector_tpetra->doExport(*localRhsVector_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    globalRhsVector_tpetra->scale(-1.);
//    globalA_JJ_xpetra->doExport(*localA_JJ_xpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::ADD);
//    globalRhsVector_tpetra.doExport(*localRhsVector_xpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
//    globalRhsVector_tpetra.scale(-1.);

//    RCP<Xpetra::CrsMatrix<double, int, int>> localA_JJ_xpetra = converter2.convertEigen2XpetraCrsMatrix(A_JJ_rowMajor, overlappingGraph_xpetra, true);
//    RCP<Xpetra::Vector<double, int, int>> localRhsVector_xpetra = converter2.convertEigen2XpetraVector(r_J, overlappingMap_xpetra);
//    globalA_JJ_xpetra->doExport(*localA_JJ_xpetra.get(), *exporter_xpetra.get(), Xpetra::CombineMode::ADD);
//    globalRhsVector_xpetra->doExport(*localRhsVector_xpetra.get(), *exporter_xpetra.get(), Xpetra::CombineMode::INSERT);
//    globalRhsVector_xpetra->scale(-1.);


    //******************************************
    //*        solve complete problem          *
    //******************************************
    Epetra_MultiVector sol = structure.solveSystem(globalA_JJ, globalRhsVector, false);
    Teuchos::RCP<Tpetra::MultiVector<double, int, int>> sol_tpetra = solveSystem_tpetra(globalA_JJ_tpetra, globalRhsVector_tpetra, false);
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

//    for (int mogi : myOwningGlobalDofIDs)
//        std::cout << rank << ": " << mogi << std::endl;

    structure.visualizeSerializedParticularSolution(sol_std, myOwningGlobalDofIDs, node2Dofs_vector, "result" + std::to_string(rank) + ".vtu", numProc);
//    structure.visualizeSolution(sol_std, "result_default" + std::to_string(rank) + ".vtu");
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



using namespace NuTo;

MeshFem QuadPatchTestMesh()
{
    /* Something like this:
     *
     *    3-----------------------2
     * /| | - _        e2       / | -->
     * /| |     -7------------6   | -->
     * /| | e4  /     e3      |   | -->
     * /| |    /              |e1 | -->
     * /| |   /    _____------5   | -->
     * /| |  4-----            \  | --> p
     * /| | /      e0           \ | -->
     * /| |/                     \| -->
     *    0-----------------------1
     *   /_\
     *   ///
     *              (c) ttitsche :)
     */
    MeshFem mesh;
    auto& n0 = mesh.Nodes.Add(Eigen::Vector2d(0, 0));
    auto& n1 = mesh.Nodes.Add(Eigen::Vector2d(10, 0));
    auto& n2 = mesh.Nodes.Add(Eigen::Vector2d(10, 10));
    auto& n3 = mesh.Nodes.Add(Eigen::Vector2d(0, 10));

    auto& n4 = mesh.Nodes.Add(Eigen::Vector2d(2, 2));
    auto& n5 = mesh.Nodes.Add(Eigen::Vector2d(8, 3));
    auto& n6 = mesh.Nodes.Add(Eigen::Vector2d(8, 7));
    auto& n7 = mesh.Nodes.Add(Eigen::Vector2d(4, 7));

    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    mesh.Elements.Add({{{n0, n1, n5, n4}, interpolation}});
    mesh.Elements.Add({{{n1, n2, n6, n5}, interpolation}});
    mesh.Elements.Add({{{n7, n6, n2, n3}, interpolation}});
    mesh.Elements.Add({{{n4, n5, n6, n7}, interpolation}});
    mesh.Elements.Add({{{n0, n4, n7, n3}, interpolation}});

    return mesh;
}


MeshFem QuadPatchTestMesh_partition0()
{
    /* Something like this:
     *
     *    3------------------------2
     * /| | - _         e2       // | -->
     * /| |     -7--------------6   | -->
     * /| | e4   /     e3      ||   | -->
     * /| |     /              ||e1 | -->
     * /| |    /         ======5    | -->
     * /| |   /     =====       \   |
     * /| |  4======             \  | --> p
     * /| | //      e0            \ | -->
     * /| |//                      \| -->
     *    0-------------------------1
     *   /_\
     *   ///
     *              (c) ttitsche :)
     */

    MeshFem mesh0;

    auto& n0 = mesh0.Nodes.Add(Eigen::Vector2d(0,0));
    auto& n2 = mesh0.Nodes.Add(Eigen::Vector2d(10,10));
    auto& n3 = mesh0.Nodes.Add(Eigen::Vector2d(0,10));
    auto& n4 = mesh0.Nodes.Add(Eigen::Vector2d(2,2));
    auto& n5 = mesh0.Nodes.Add(Eigen::Vector2d(8,3));
    auto& n6 = mesh0.Nodes.Add(Eigen::Vector2d(8,7));
    auto& n7 = mesh0.Nodes.Add(Eigen::Vector2d(4,7));

    const auto& interpolation0 = mesh0.CreateInterpolation(InterpolationQuadLinear(2));

    mesh0.Elements.Add({{{n7, n6, n2, n3}, interpolation0}});
    mesh0.Elements.Add({{{n4, n5, n6, n7}, interpolation0}});
    mesh0.Elements.Add({{{n0, n4, n7, n3}, interpolation0}});

    return mesh0;
}


MeshFem QuadPatchTestMesh_partition1()
{
    /* Something like this:
     *
     *    3------------------------2
     * /| | - _         e2       // | -->
     * /| |     -7--------------6   | -->
     * /| | e4   /     e3      ||   | -->
     * /| |     /              ||e1 | -->
     * /| |    /         ======5    | -->
     * /| |   /     =====       \   |
     * /| |  4======             \  | --> p
     * /| | //      e0            \ | -->
     * /| |//                      \| -->
     *    0-------------------------1
     *   /_\
     *   ///
     *              (c) ttitsche :)
     */

    MeshFem mesh1;

    auto& m0 = mesh1.Nodes.Add(Eigen::Vector2d(0,0));
    auto& m1 = mesh1.Nodes.Add(Eigen::Vector2d(10,0));
    auto& m2 = mesh1.Nodes.Add(Eigen::Vector2d(10,10));
    auto& m4 = mesh1.Nodes.Add(Eigen::Vector2d(2,2));
    auto& m5 = mesh1.Nodes.Add(Eigen::Vector2d(8,3));
    auto& m6 = mesh1.Nodes.Add(Eigen::Vector2d(8,7));

    const auto& interpolation1 = mesh1.CreateInterpolation(InterpolationQuadLinear(2));

    mesh1.Elements.Add({{{m0, m1, m5, m4}, interpolation1}});
    mesh1.Elements.Add({{{m1, m2, m6, m5}, interpolation1}});

    return mesh1;
}


struct DofInfo
{
    DofContainer<int> numIndependentDofs;
    DofContainer<int> numDependentDofs;

    std::vector<int> activeDofs;
    std::vector<int> dependentDofs;
};

DofInfo ManualDofNumbering(MeshFem* rMesh, DofType dof)
{
    // some manual dof numbering ...
    NuTo::Groups::Group<NodeSimple> allNodes = rMesh->NodesTotal(dof);
    NuTo::Groups::Group<NodeSimple> nodesConstrainedInX = rMesh->NodesAtAxis(eDirection::X, dof);
    NuTo::Groups::Group<NodeSimple> nodesConstrainedInY = NuTo::Groups::Group<NodeSimple>(rMesh->NodeAtCoordinate(Eigen::Vector2d(0, 0), dof));

    NuTo::Groups::Group<NodeSimple> nodesUnconstrainedInX = Difference(allNodes, nodesConstrainedInX);
    NuTo::Groups::Group<NodeSimple> nodesUnconstrainedInY = Difference(allNodes, nodesConstrainedInY);

    DofInfo dofInfo;
    int dofNumber = 0;

    for (auto& node : nodesUnconstrainedInX)
    {
        node.SetDofNumber(0, dofNumber);
        dofInfo.activeDofs.push_back(dofNumber);
        dofNumber++;
    }
    for (auto& node : nodesUnconstrainedInY)
    {
        node.SetDofNumber(1, dofNumber);
        dofInfo.activeDofs.push_back(dofNumber);
        dofNumber++;
    }

    for (auto& node : nodesConstrainedInX)
    {
        node.SetDofNumber(0, dofNumber);
        dofInfo.dependentDofs.push_back(dofNumber);
        dofNumber++;
    }
    for (auto& node : nodesConstrainedInY)
    {
        node.SetDofNumber(1, dofNumber);
        dofInfo.dependentDofs.push_back(dofNumber);
        dofNumber++;
    }


    dofInfo.numIndependentDofs[dof] = nodesUnconstrainedInX.Size() + nodesUnconstrainedInY.Size();
    dofInfo.numDependentDofs[dof] = nodesConstrainedInX.Size() + nodesConstrainedInY.Size();
    return dofInfo;
}

std::vector<DofInfo> dofNumberingSubMeshes(MeshFem* rMesh, MeshFem* rSubMesh0, MeshFem* rSubMesh1, DofType rDof)
{
    std::vector<DofInfo> numberings(2);

    numberings[0] = ManualDofNumbering(rSubMesh0, rDof);
    numberings[1] = ManualDofNumbering(rSubMesh1, rDof);

    return numberings;
}


std::vector<std::map<int, int>> local2GlobalNumbering(MeshFem* rGlobalMesh, MeshFem* rLocalMesh0, MeshFem* rLocalMesh1, DofType rDof)
{
    std::vector<std::map<int, int>> local2Global(2);

    NodeSimple n0 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(0,0), rDof);
    NodeSimple& n00 = rLocalMesh0->NodeAtCoordinate(Eigen::Vector2d(0,0), rDof);
    NodeSimple& n10 = rLocalMesh1->NodeAtCoordinate(Eigen::Vector2d(0,0), rDof);
    local2Global[0][n00.GetDofNumber(0)] = n0.GetDofNumber(0);
    local2Global[0][n00.GetDofNumber(1)] = n0.GetDofNumber(1);
    local2Global[1][n10.GetDofNumber(0)] = n0.GetDofNumber(0);
    local2Global[1][n10.GetDofNumber(1)] = n0.GetDofNumber(1);

    NodeSimple n1 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(10,0), rDof);
    NodeSimple& n11 = rLocalMesh1->NodeAtCoordinate(Eigen::Vector2d(10,0), rDof);
    local2Global[1][n11.GetDofNumber(0)] = n1.GetDofNumber(0);
    local2Global[1][n11.GetDofNumber(1)] = n1.GetDofNumber(1);

    NodeSimple n2 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(10,10), rDof);
    NodeSimple& n02 = rLocalMesh0->NodeAtCoordinate(Eigen::Vector2d(10,10), rDof);
    NodeSimple& n12 = rLocalMesh1->NodeAtCoordinate(Eigen::Vector2d(10,10), rDof);
    local2Global[0][n02.GetDofNumber(0)] = n2.GetDofNumber(0);
    local2Global[0][n02.GetDofNumber(1)] = n2.GetDofNumber(1);
    local2Global[1][n12.GetDofNumber(0)] = n2.GetDofNumber(0);
    local2Global[1][n12.GetDofNumber(1)] = n2.GetDofNumber(1);

    NodeSimple n3 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(0,10), rDof);
    NodeSimple& n03 = rLocalMesh0->NodeAtCoordinate(Eigen::Vector2d(0,10), rDof);
    local2Global[0][n03.GetDofNumber(0)] = n3.GetDofNumber(0);
    local2Global[0][n03.GetDofNumber(1)] = n3.GetDofNumber(1);

    NodeSimple n4 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(2,2), rDof);
    NodeSimple& n04 = rLocalMesh0->NodeAtCoordinate(Eigen::Vector2d(2,2), rDof);
    NodeSimple& n14 = rLocalMesh1->NodeAtCoordinate(Eigen::Vector2d(2,2), rDof);
    local2Global[0][n04.GetDofNumber(0)] = n4.GetDofNumber(0);
    local2Global[0][n04.GetDofNumber(1)] = n4.GetDofNumber(1);
    local2Global[1][n14.GetDofNumber(0)] = n4.GetDofNumber(0);
    local2Global[1][n14.GetDofNumber(1)] = n4.GetDofNumber(1);

    NodeSimple n5 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(8,3), rDof);
    NodeSimple& n05 = rLocalMesh0->NodeAtCoordinate(Eigen::Vector2d(8,3), rDof);
    NodeSimple& n15 = rLocalMesh1->NodeAtCoordinate(Eigen::Vector2d(8,3), rDof);
    local2Global[0][n05.GetDofNumber(0)] = n5.GetDofNumber(0);
    local2Global[0][n05.GetDofNumber(1)] = n5.GetDofNumber(1);
    local2Global[1][n15.GetDofNumber(0)] = n5.GetDofNumber(0);
    local2Global[1][n15.GetDofNumber(1)] = n5.GetDofNumber(1);

    NodeSimple n6 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(8,7), rDof);
    NodeSimple& n06 = rLocalMesh0->NodeAtCoordinate(Eigen::Vector2d(8,7), rDof);
    NodeSimple& n16 = rLocalMesh1->NodeAtCoordinate(Eigen::Vector2d(8,7), rDof);
    local2Global[0][n06.GetDofNumber(0)] = n6.GetDofNumber(0);
    local2Global[0][n06.GetDofNumber(1)] = n6.GetDofNumber(1);
    local2Global[1][n16.GetDofNumber(0)] = n6.GetDofNumber(0);
    local2Global[1][n16.GetDofNumber(1)] = n6.GetDofNumber(1);

    NodeSimple n7 = rGlobalMesh->NodeAtCoordinate(Eigen::Vector2d(4,7), rDof);
    NodeSimple& n07 = rLocalMesh0->NodeAtCoordinate(Eigen::Vector2d(4,7), rDof);
    local2Global[0][n07.GetDofNumber(0)] = n7.GetDofNumber(0);
    local2Global[0][n07.GetDofNumber(1)] = n7.GetDofNumber(1);

    return local2Global;
}

std::vector<int> getAllDofNumbers(std::vector<std::map<int, int>> rLocal2GlobalDofs, int rRank)
{
    StructureMesh structure(2);
    std::vector<int> dofNumbers = structure.map2Vector(rLocal2GlobalDofs[rRank]);

    return dofNumbers;
}

std::vector<int> getOwningDofNumbers(std::vector<std::map<int, int>> rLocal2GlobalDofs, int rRank)
{
    std::vector<int> dofNumbers;
    if (rRank == 0)
    {
        StructureMesh structure(2);
        dofNumbers = structure.map2Vector(rLocal2GlobalDofs[rRank]);
    }
    else
    {
        dofNumbers.push_back(rLocal2GlobalDofs[rRank][0]);
        dofNumbers.push_back(rLocal2GlobalDofs[rRank][5]);
    }

    return dofNumbers;
}

std::vector<int> getGlobalOwningActiveDofNumbers(std::vector<std::map<int,int>> rLocal2GlobalDofs, int rRank, std::vector<int> localActiveDofs)
{
    std::vector<int> globalActiveDofs;

    if (rRank == 0)
    {
        for (int localDof : localActiveDofs)
        {
            globalActiveDofs.push_back(rLocal2GlobalDofs[rRank][localDof]);
        }
    }
    else
    {
        globalActiveDofs.push_back(rLocal2GlobalDofs[rRank][0]);
        globalActiveDofs.push_back(rLocal2GlobalDofs[rRank][5]);
    }
    return globalActiveDofs;
}


void run_simpleAssembler_test()
{

    MeshFem mesh = QuadPatchTestMesh();
    DofType displ("displacements", 2);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    AddDofInterpolation(&mesh, displ, interpolation);
    DofInfo dofInfo = ManualDofNumbering(&mesh, displ);


    // ************************************************************************
    //                 add continuum cells
    // ************************************************************************
    constexpr double E = 20000;
    constexpr double nu = 0.2;
    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::TimeDependent::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);

    boost::ptr_vector<CellInterface> cellContainer;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

    NuTo::Groups::Group<CellInterface> cellGroup;
    for (auto& element : mesh.Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, momentumBalance));
        cellGroup.Add(cellContainer.back());
    }

    // ************************************************************************
    //                 add boundary cells
    // ************************************************************************

    // manually add the boundary element
    const auto& interpolationBc = mesh.CreateInterpolation(InterpolationTrussLinear(2));

    // extract existing nodes
    auto boundaryCoordNodes = mesh.NodesAtAxis(eDirection::X, 10);
    NodeSimple& nc1 = *boundaryCoordNodes.begin();
    NodeSimple& nc2 = *(boundaryCoordNodes.begin() + 1);

    auto boundaryDisplNodes = mesh.NodesAtAxis(eDirection::X, displ, 10);
    NodeSimple& nd1 = *boundaryDisplNodes.begin();
    NodeSimple& nd2 = *(boundaryDisplNodes.begin() + 1);

    // add the boundary element
    auto& boundaryElement = mesh.Elements.Add({{{nc1, nc2}, interpolationBc}});
    boundaryElement.AddDofElement(displ, {{nd1, nd2}, interpolationBc});

    IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBC(1, 0);
    Integrands::TimeDependent::NeumannBc<2> neumannBc(displ, pressureBC);

    cellContainer.push_back(new Cell(boundaryElement, integrationTypeBc, neumannBc));
    cellGroup.Add(cellContainer.back());

    // ************************************************************************
    //                  assemble and solve
    // ************************************************************************
    SimpleAssembler assembler(dofInfo.numIndependentDofs, dofInfo.numDependentDofs);

    auto gradient = assembler.BuildVector(cellGroup, {&displ}, Integrands::TimeDependent::Gradient());
    auto hessian = assembler.BuildMatrix(cellGroup, {&displ}, Integrands::TimeDependent::Hessian0());

    Eigen::MatrixXd hessianDense(hessian.JJ(displ, displ));
    Eigen::VectorXd newDisplacements = hessianDense.ldlt().solve(gradient.J[displ]);

    // merge dof values
    int numUnconstrainedDofs = dofInfo.numIndependentDofs[displ];
    for (auto& node : mesh.NodesTotal(displ))
    {
        int dofX = node.GetDofNumber(0);
        int dofY = node.GetDofNumber(1);

        if (dofX < numUnconstrainedDofs)
            node.SetValue(0, newDisplacements[dofX]);
        if (dofY < numUnconstrainedDofs)
            node.SetValue(1, newDisplacements[dofY]);
    }

    for (auto& node : mesh.NodesTotal(displ))
    {
        Eigen::VectorXd nodValues = node.GetValues();
        for (int i = 0; i < nodValues.rows(); ++i)
        {
            std::cout << nodValues(i) << "  ";
        }
        std::cout << std::endl;
    }

    ConversionTools converter2;
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> A_JJ_tpetra = converter2.convertEigen2TpetraCrsMatrix(hessian.JJ(displ, displ), true, true);
    Teuchos::RCP<Tpetra::Vector<double, int, int>> rhsVector_tpetra = converter2.convertEigen2TpetraVector(gradient.J[displ], true);

    if (converter2.getComm_tpetra()->getRank() == 0)
    {
        std::cout << hessian.JJ(displ, displ).toDense() << std::endl;
        std::cout << gradient.J[displ] << std::endl;
    }

    Teuchos::RCP<Tpetra::Vector<double, int, int>> rhsV = rcp(new Tpetra::Vector<double, int, int>(A_JJ_tpetra->getRangeMap()));
    Teuchos::RCP<Tpetra::MultiVector<double, int, int>> sol_tpetra = solveSystem_tpetra(A_JJ_tpetra, rhsVector_tpetra, false);

}



NuTo::Groups::Group<CellInterface> setUpMeshCells(MeshFem* rMesh, DofType rDofType, boost::ptr_vector<CellInterface>* cellContainer, NuTo::Groups::Group<CellInterface>* cellGroup)
{
    // ************************************************************************
    //                 add continuum cells
    // ************************************************************************
    constexpr double E = 20000;
    constexpr double nu = 0.2;
    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::TimeDependent::MomentumBalance<2> momentumBalance(rDofType, linearElasticLaw);

//    boost::ptr_vector<CellInterface> cellContainer;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

//    NuTo::Groups::Group<CellInterface> cellGroup;
    for (auto& element : rMesh->Elements)
    {
        cellContainer->push_back(new Cell(element, integrationType, momentumBalance));
        cellGroup->Add(cellContainer->back());
    }

    // ************************************************************************
    //                 add boundary cells
    // ************************************************************************

    // manually add the boundary element
    const auto& interpolationBc = rMesh->CreateInterpolation(InterpolationTrussLinear(2));

    // extract existing nodes
    NuTo::Groups::Group<NodeSimple> boundaryCoordNodes = rMesh->NodesAtAxis(eDirection::X, 10);
    if (boundaryCoordNodes.Size() == 2)
    {
        NodeSimple& nc1 = *boundaryCoordNodes.begin();
        NodeSimple& nc2 = *(boundaryCoordNodes.begin() + 1);

        auto boundaryDisplNodes = rMesh->NodesAtAxis(eDirection::X, rDofType, 10);
        NodeSimple& nd1 = *boundaryDisplNodes.begin();
        NodeSimple& nd2 = *(boundaryDisplNodes.begin() + 1);

        // add the boundary element
        auto& boundaryElement = rMesh->Elements.Add({{{nc1, nc2}, interpolationBc}});
        boundaryElement.AddDofElement(rDofType, {{nd1, nd2}, interpolationBc});

        IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
        Eigen::Vector2d pressureBC(1, 0);
        Integrands::TimeDependent::NeumannBc<2> neumannBc(rDofType, pressureBC);

        cellContainer->push_back(new Cell(boundaryElement, integrationTypeBc, neumannBc));
        cellGroup->Add(cellContainer->back());
    }

    return *cellGroup;
}


void run_Assembler_test()
{
    RCP<const Teuchos::Comm<int>> commTeuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    int rank = commTeuchos->getRank();
    MeshFem mesh = QuadPatchTestMesh();
    DofType displ("displacements", 2);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    AddDofInterpolation(&mesh, displ, interpolation);
    DofInfo dofInfo = ManualDofNumbering(&mesh, displ);

    //-------------- mesh partitioning
    std::vector<MeshFem> subMeshes(2);
    MeshFem subMesh0 = QuadPatchTestMesh_partition0();
    MeshFem subMesh1 = QuadPatchTestMesh_partition1();

    const auto& interpolation0 = subMesh0.CreateInterpolation(InterpolationQuadLinear(2));
    AddDofInterpolation(&subMesh0, displ, interpolation0);
    const auto& interpolation1 = subMesh1.CreateInterpolation(InterpolationQuadLinear(2));
    AddDofInterpolation(&subMesh1, displ, interpolation1);

    std::vector<DofInfo> subMeshNumbering = dofNumberingSubMeshes(&mesh, &subMesh0, &subMesh1, displ);
    std::vector<std::map<int, int>> local2GlobalDofs = local2GlobalNumbering(&mesh, &subMesh0, &subMesh1, displ);
    subMeshes[0] = std::move(subMesh0);
    subMeshes[1] = std::move(subMesh1);

    //--------- do the following steps for every part of the mesh

    DofInfo currDofInfo = subMeshNumbering[rank];

    // ************************************************************************
    //                 add continuum cells
    // ************************************************************************
    constexpr double E = 20000;
    constexpr double nu = 0.2;
    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::TimeDependent::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);

    boost::ptr_vector<CellInterface> cellContainer;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

    NuTo::Groups::Group<CellInterface> cellGroup;
    for (auto& element : subMeshes[rank].Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, momentumBalance));
        cellGroup.Add(cellContainer.back());
    }

    // ************************************************************************
    //                 add boundary cells
    // ************************************************************************

    // manually add the boundary element
    const auto& interpolationBc = subMeshes[rank].CreateInterpolation(InterpolationTrussLinear(2));

    IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBC(1, 0);
    Integrands::TimeDependent::NeumannBc<2> neumannBc(displ, pressureBC);

    // extract existing nodes
    NuTo::Groups::Group<NodeSimple> boundaryCoordNodes = subMeshes[rank].NodesAtAxis(eDirection::X, 10);
    if (boundaryCoordNodes.Size() == 2)
    {
        NodeSimple& nc1 = *boundaryCoordNodes.begin();
        NodeSimple& nc2 = *(boundaryCoordNodes.begin() + 1);

        auto boundaryDisplNodes = subMeshes[rank].NodesAtAxis(eDirection::X, displ, 10);
        NodeSimple& nd1 = *boundaryDisplNodes.begin();
        NodeSimple& nd2 = *(boundaryDisplNodes.begin() + 1);

        // add the boundary element
        auto& boundaryElement = subMeshes[rank].Elements.Add({{{nc1, nc2}, interpolationBc}});
        boundaryElement.AddDofElement(displ, {{nd1, nd2}, interpolationBc});

        cellContainer.push_back(new Cell(boundaryElement, integrationTypeBc, neumannBc));
        cellGroup.Add(cellContainer.back());
    }

    // ************************************************************************
    //                  assemble
    // ************************************************************************
    SimpleAssembler assembler(currDofInfo.numIndependentDofs, currDofInfo.numDependentDofs);

    auto gradient = assembler.BuildVector(cellGroup, {&displ}, Integrands::TimeDependent::Gradient());
    auto hessian = assembler.BuildMatrix(cellGroup, {&displ}, Integrands::TimeDependent::Hessian0());
    Eigen::SparseMatrix<double> A_JJ = hessian.JJ(displ,displ);
    Eigen::VectorXd r_J = gradient.J[displ];

    //-----------

    //******************************************
    //*     create overlapping index map       *
    //******************************************
    std::vector<int> myGlobalDofIDs = getAllDofNumbers(local2GlobalDofs, rank);
    int* myGlobalDofIDs_arr = &myGlobalDofIDs[0];
    RCP<Tpetra::Map<int, int>> overlappingMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myGlobalDofIDs_arr, myGlobalDofIDs.size(), 0, commTeuchos));


    //******************************************
    //*       create owning index map          *
    //******************************************
    std::vector<int> myOwningGlobalActiveDofIDs = getGlobalOwningActiveDofNumbers(local2GlobalDofs, rank, currDofInfo.activeDofs);
    int* myOwningGlobalActiveDofIDs_arr = &myOwningGlobalActiveDofIDs[0];
    RCP<Tpetra::Map<int, int>> owningMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myOwningGlobalActiveDofIDs_arr, myOwningGlobalActiveDofIDs.size(), 0, commTeuchos));

//    std::ostream &out = std::cout;
//    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
//    *fos << "OwningMap :" << std::endl;
//    owningMap_tpetra->describe(*fos,Teuchos::VERB_EXTREME);
//    *fos << std::endl;
//    *fos << "OverlappingMap :" << std::endl;
//    overlappingMap_tpetra->describe(*fos,Teuchos::VERB_EXTREME);
//    *fos << std::endl;

    //******************************************
    //*         create index graphs            *
    //******************************************
    int maxNonZeros = 18;   //got by interpolation type (order), e.g. QUAD2 + EQUIDISTANT1 => 18
    RCP<Tpetra::CrsGraph<int, int>> owningGraph_tpetra = rcp(new Tpetra::CrsGraph<int,int>(owningMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
    RCP<Tpetra::CrsGraph<int, int>> overlappingGraph_tpetra = rcp(new Tpetra::CrsGraph<int, int>(overlappingMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
    std::vector<int> columnIndices;
    for (int k=0; k<A_JJ.outerSize(); ++k)
    {
        columnIndices.clear();
        for (Eigen::SparseMatrix<double>::InnerIterator it(A_JJ,k); it; ++it)
        {
            // describe position of entries
            overlappingGraph_tpetra->insertGlobalIndices(myGlobalDofIDs_arr[it.row()], 1, &myGlobalDofIDs_arr[it.col()]);
        }
    }

    //******************************************
    //*   define inter-process communication   *
    //*      for local-to-global indices       *
    //******************************************
    RCP<const Tpetra::Export<int, int>> exporter_tpetra = rcp(new Tpetra::Export<int, int>(overlappingMap_tpetra, owningMap_tpetra));
    overlappingGraph_tpetra->fillComplete();
    owningGraph_tpetra->doExport(*overlappingGraph_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    owningGraph_tpetra->fillComplete();


    //******************************************
    //*  initialize Trilinos matrix and vector *
    //******************************************
    RCP<Tpetra::CrsMatrix<double, int, int>> globalA_JJ_tpetra = rcp(new Tpetra::CrsMatrix<double, int, int>(owningGraph_tpetra));
    RCP<Tpetra::Vector<double, int, int>> globalRhsVector_tpetra = rcp(new Tpetra::Vector<double, int, int>(owningMap_tpetra));
    globalRhsVector_tpetra->putScalar(0.0);

    Eigen::SparseMatrix<double, Eigen::RowMajor> A_JJ_rowMajor(A_JJ);

    //******************************************
    //*    conversion from NuTo to Trilinos    *
    //******************************************
    ConversionTools converter2;
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> localA_JJ_tpetra = converter2.convertEigen2TpetraCrsMatrix(A_JJ_rowMajor, overlappingGraph_tpetra, true);
    Teuchos::RCP<Tpetra::Vector<double, int, int>> localRhsVector_tpetra = converter2.convertEigen2TpetraVector(r_J, overlappingMap_tpetra);
    globalA_JJ_tpetra->doExport(*localA_JJ_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::ADD);
    globalRhsVector_tpetra->doExport(*localRhsVector_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    globalRhsVector_tpetra->scale(-1.);

    //******************************************
    //*        solve complete problem          *
    //******************************************
    Teuchos::RCP<Tpetra::MultiVector<double, int, int>> sol_tpetra = solveSystem_tpetra(globalA_JJ_tpetra, globalRhsVector_tpetra, false);
}

std::vector<ZoltanMesh> partitionMesh(int argc, char** argv, ZoltanMesh* rMesh, int rRank, int rNumProcs = 2)
{
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
    int *importProcs, *importToPart, *exportProcs, *exportToPart;

    struct Zoltan_Struct *zz;
    int zoltanInfo = 0;
    float version = 0;
    zoltanInfo = Zoltan_Initialize(argc, argv, &version);


    if (zoltanInfo != ZOLTAN_OK){
        std::cout << "ERROR: Something went wrong with Zoltan" << std::endl;
//        MPI_Finalize();
        exit(0);
    }

    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */
    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
    Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); /* use Zoltan default vertex weights */
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */

    /* PHG parameters  - see the Zoltan User's Guide for many more
    *   (The "REPARTITION" approach asks Zoltan to create a partitioning that is
    *    better but is not too far from the current partitioning, rather than partitioning
    *    from scratch.  It may be faster but of lower quality that LB_APPROACH=PARTITION.)
    */

//    Zoltan_Set_Param(zz, "LB_APPROACH", "REPARTITION");
    Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");

    /* Application defined query functions */
    Zoltan_Set_Num_Obj_Fn(zz, rMesh->get_number_of_localElements, rMesh);
    Zoltan_Set_Obj_List_Fn(zz, rMesh->get_localElement_list, rMesh);
    Zoltan_Set_HG_Size_CS_Fn(zz, rMesh->get_number_of_localNodes_localPins, rMesh);
    Zoltan_Set_HG_CS_Fn(zz, rMesh->get_hypergraph, rMesh);

    zoltanInfo = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
                             &changes,        /* 1 if partitioning was changed, 0 otherwise */
                             &numGidEntries,  /* Number of integers used for a global ID */
                             &numLidEntries,  /* Number of integers used for a local ID */
                             &numImport,      /* Number of vertices to be sent to me */
                             &importGlobalGids,  /* Global IDs of vertices to be sent to me */
                             &importLocalGids,   /* Local IDs of vertices to be sent to me */
                             &importProcs,    /* Process rank for source of each incoming vertex */
                             &importToPart,   /* New partition for each incoming vertex */
                             &numExport,      /* Number of vertices I must send to other processes*/
                             &exportGlobalGids,  /* Global IDs of the vertices I must send */
                             &exportLocalGids,   /* Local IDs of the vertices I must send */
                             &exportProcs,    /* Process to which I send each of the vertices */
                             &exportToPart);  /* Partition to which each vertex will belong */

    if (zoltanInfo != ZOLTAN_OK){
        std::cout << "ERROR: Something went wrong with Zoltan" << std::endl;
//        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

//    if (myRank == 0)
//    {
        std::cout << "Partition summary\n-----------------" << std::endl;
        std::cout << "changes: " << changes << std::endl;
        std::cout << "numGidEntries: " << numGidEntries << std::endl;
        std::cout << "numLidEntries: " << numLidEntries << std::endl;
        std::cout << "numImport: " << numImport << std::endl;
        std::cout << "importGlobalGids: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importGlobalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importLocalGids: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importLocalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importProcs: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importProcs[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "importToPart: " << std::endl;
        for (int i = 0; i < numImport; ++i)
        {
            std::cout << importToPart[i] << ", ";
        }
        std::cout << "numExport: " << numExport << std::endl;
        std::cout << "exportGlobalGids: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportGlobalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportLocalGids: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportLocalGids[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportProcs: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportProcs[i] << ", ";
        }
        std::cout << std::endl;
        std::cout << "exportToPart: " << std::endl;
        for (int i = 0; i < numExport; ++i)
        {
            std::cout << exportToPart[i] << ", ";
        }
        std::cout << std::endl;
//    }

    int parts[rMesh->numLocalElements] = {rRank};

    for (int i=0; i < numExport; ++i){
        parts[exportLocalGids[i]] = exportToPart[i];
    }

    int partAssign[rMesh->numGlobalElements], allPartAssign[rMesh->numGlobalElements];

    memset(partAssign, 0, sizeof(int) * rMesh->numGlobalElements);

    for (int i=0; i < rMesh->numLocalElements; ++i){
        partAssign[rMesh->myGlobalElementIDs[i]] = parts[i];
    }
    MPI_Reduce(partAssign, allPartAssign, rMesh->numGlobalElements, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    std::cout << "Local elements of Proc [" << rRank << "]: ";
    for (int part : partAssign)
    {
        std::cout << part << " ";
    }
    std::cout << std::endl;

    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,
                        &exportProcs, &exportToPart);

    Zoltan_Destroy(&zz);
}

void run_Assembler_Zoltan_test(int argc, char** argv)
{
    RCP<const Teuchos::Comm<int>> commTeuchos = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

    int rank = commTeuchos->getRank();
    ZoltanMesh mesh = ZoltanMesh::create_2D_mesh();
    DofType displ("displacements", 2);
    const auto& interpolation = mesh.CreateInterpolation(InterpolationQuadLinear(2));

    AddDofInterpolation(&mesh, displ, interpolation);
    DofInfo dofInfo = ManualDofNumbering(&mesh, displ);

    //-------------- mesh partitioning
    std::vector<MeshFem> subMeshes(2);
    MeshFem subMesh0 = QuadPatchTestMesh_partition0();
    MeshFem subMesh1 = QuadPatchTestMesh_partition1();

    std::cout << "Rank[" << commTeuchos->getRank() << "] : localElements (before) = " << mesh.numLocalElements << std::endl;
    partitionMesh(argc, argv ,&mesh, rank);
    std::cout << "Rank[" << commTeuchos->getRank() << "] : localElements (after) = " << mesh.numLocalElements << std::endl;

    const auto& interpolation0 = subMesh0.CreateInterpolation(InterpolationQuadLinear(2));
    AddDofInterpolation(&subMesh0, displ, interpolation0);
    const auto& interpolation1 = subMesh1.CreateInterpolation(InterpolationQuadLinear(2));
    AddDofInterpolation(&subMesh1, displ, interpolation1);

    std::vector<DofInfo> subMeshNumbering = dofNumberingSubMeshes(&mesh, &subMesh0, &subMesh1, displ);
    std::vector<std::map<int, int>> local2GlobalDofs = local2GlobalNumbering(&mesh, &subMesh0, &subMesh1, displ);
    subMeshes[0] = std::move(subMesh0);
    subMeshes[1] = std::move(subMesh1);

    //--------- do the following steps for every part of the mesh

    DofInfo currDofInfo = subMeshNumbering[rank];

    // ************************************************************************
    //                 add continuum cells
    // ************************************************************************
    constexpr double E = 20000;
    constexpr double nu = 0.2;
    Laws::LinearElastic<2> linearElasticLaw(E, nu, ePlaneState::PLANE_STRESS);
    Integrands::TimeDependent::MomentumBalance<2> momentumBalance(displ, linearElasticLaw);

    boost::ptr_vector<CellInterface> cellContainer;
    IntegrationTypeTensorProduct<2> integrationType(2, eIntegrationMethod::GAUSS);

    NuTo::Groups::Group<CellInterface> cellGroup;
    for (auto& element : subMeshes[rank].Elements)
    {
        cellContainer.push_back(new Cell(element, integrationType, momentumBalance));
        cellGroup.Add(cellContainer.back());
    }

    // ************************************************************************
    //                 add boundary cells
    // ************************************************************************

    // manually add the boundary element
    const auto& interpolationBc = subMeshes[rank].CreateInterpolation(InterpolationTrussLinear(2));

    IntegrationTypeTensorProduct<1> integrationTypeBc(1, eIntegrationMethod::GAUSS);
    Eigen::Vector2d pressureBC(1, 0);
    Integrands::TimeDependent::NeumannBc<2> neumannBc(displ, pressureBC);

    // extract existing nodes
    NuTo::Groups::Group<NodeSimple> boundaryCoordNodes = subMeshes[rank].NodesAtAxis(eDirection::X, 10);
    if (boundaryCoordNodes.Size() == 2)
    {
        NodeSimple& nc1 = *boundaryCoordNodes.begin();
        NodeSimple& nc2 = *(boundaryCoordNodes.begin() + 1);

        auto boundaryDisplNodes = subMeshes[rank].NodesAtAxis(eDirection::X, displ, 10);
        NodeSimple& nd1 = *boundaryDisplNodes.begin();
        NodeSimple& nd2 = *(boundaryDisplNodes.begin() + 1);

        // add the boundary element
        auto& boundaryElement = subMeshes[rank].Elements.Add({{{nc1, nc2}, interpolationBc}});
        boundaryElement.AddDofElement(displ, {{nd1, nd2}, interpolationBc});

        cellContainer.push_back(new Cell(boundaryElement, integrationTypeBc, neumannBc));
        cellGroup.Add(cellContainer.back());
    }

    // ************************************************************************
    //                  assemble
    // ************************************************************************
    SimpleAssembler assembler(currDofInfo.numIndependentDofs, currDofInfo.numDependentDofs);

    auto gradient = assembler.BuildVector(cellGroup, {&displ}, Integrands::TimeDependent::Gradient());
    auto hessian = assembler.BuildMatrix(cellGroup, {&displ}, Integrands::TimeDependent::Hessian0());
    Eigen::SparseMatrix<double> A_JJ = hessian.JJ(displ,displ);
    Eigen::VectorXd r_J = gradient.J[displ];

    //-----------

    //******************************************
    //*     create overlapping index map       *
    //******************************************
    std::vector<int> myGlobalDofIDs = getAllDofNumbers(local2GlobalDofs, rank);
    int* myGlobalDofIDs_arr = &myGlobalDofIDs[0];
    RCP<Tpetra::Map<int, int>> overlappingMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myGlobalDofIDs_arr, myGlobalDofIDs.size(), 0, commTeuchos));


    //******************************************
    //*       create owning index map          *
    //******************************************
    std::vector<int> myOwningGlobalActiveDofIDs = getGlobalOwningActiveDofNumbers(local2GlobalDofs, rank, currDofInfo.activeDofs);
    int* myOwningGlobalActiveDofIDs_arr = &myOwningGlobalActiveDofIDs[0];
    RCP<Tpetra::Map<int, int>> owningMap_tpetra = rcp(new Tpetra::Map<int, int>(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), myOwningGlobalActiveDofIDs_arr, myOwningGlobalActiveDofIDs.size(), 0, commTeuchos));

//    std::ostream &out = std::cout;
//    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
//    *fos << "OwningMap :" << std::endl;
//    owningMap_tpetra->describe(*fos,Teuchos::VERB_EXTREME);
//    *fos << std::endl;
//    *fos << "OverlappingMap :" << std::endl;
//    overlappingMap_tpetra->describe(*fos,Teuchos::VERB_EXTREME);
//    *fos << std::endl;

    //******************************************
    //*         create index graphs            *
    //******************************************
    int maxNonZeros = 18;   //got by interpolation type (order), e.g. QUAD2 + EQUIDISTANT1 => 18
    RCP<Tpetra::CrsGraph<int, int>> owningGraph_tpetra = rcp(new Tpetra::CrsGraph<int,int>(owningMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
    RCP<Tpetra::CrsGraph<int, int>> overlappingGraph_tpetra = rcp(new Tpetra::CrsGraph<int, int>(overlappingMap_tpetra, maxNonZeros, Tpetra::ProfileType::DynamicProfile));
    std::vector<int> columnIndices;
    for (int k=0; k<A_JJ.outerSize(); ++k)
    {
        columnIndices.clear();
        for (Eigen::SparseMatrix<double>::InnerIterator it(A_JJ,k); it; ++it)
        {
            // describe position of entries
            overlappingGraph_tpetra->insertGlobalIndices(myGlobalDofIDs_arr[it.row()], 1, &myGlobalDofIDs_arr[it.col()]);
        }
    }

    //******************************************
    //*   define inter-process communication   *
    //*      for local-to-global indices       *
    //******************************************
    RCP<const Tpetra::Export<int, int>> exporter_tpetra = rcp(new Tpetra::Export<int, int>(overlappingMap_tpetra, owningMap_tpetra));
    overlappingGraph_tpetra->fillComplete();
    owningGraph_tpetra->doExport(*overlappingGraph_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    owningGraph_tpetra->fillComplete();


    //******************************************
    //*  initialize Trilinos matrix and vector *
    //******************************************
    RCP<Tpetra::CrsMatrix<double, int, int>> globalA_JJ_tpetra = rcp(new Tpetra::CrsMatrix<double, int, int>(owningGraph_tpetra));
    RCP<Tpetra::Vector<double, int, int>> globalRhsVector_tpetra = rcp(new Tpetra::Vector<double, int, int>(owningMap_tpetra));
    globalRhsVector_tpetra->putScalar(0.0);

    Eigen::SparseMatrix<double, Eigen::RowMajor> A_JJ_rowMajor(A_JJ);

    //******************************************
    //*    conversion from NuTo to Trilinos    *
    //******************************************
    ConversionTools converter2;
    Teuchos::RCP<Tpetra::CrsMatrix<double, int, int>> localA_JJ_tpetra = converter2.convertEigen2TpetraCrsMatrix(A_JJ_rowMajor, overlappingGraph_tpetra, true);
    Teuchos::RCP<Tpetra::Vector<double, int, int>> localRhsVector_tpetra = converter2.convertEigen2TpetraVector(r_J, overlappingMap_tpetra);
    globalA_JJ_tpetra->doExport(*localA_JJ_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::ADD);
    globalRhsVector_tpetra->doExport(*localRhsVector_tpetra.get(), *exporter_tpetra.get(), Tpetra::CombineMode::INSERT);
    globalRhsVector_tpetra->scale(-1.);

    //******************************************
    //*        solve complete problem          *
    //******************************************
    Teuchos::RCP<Tpetra::MultiVector<double, int, int>> sol_tpetra = solveSystem_tpetra(globalA_JJ_tpetra, globalRhsVector_tpetra, false);
}


int main(int argc, char **argv)
{
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

//    Epetra_MpiComm comm(MPI_COMM_WORLD);
//    run_generation_test();
    std::string fileNames = "../meshes/origMesh_200.mff";
//    fileNames = "../meshes/jsonMesh_Test_NodeIDs.mesh";

//    run_mesh_test(comm, fileNames);
//    run_simpleAssembler_test();
//    run_Assembler_test();
    run_Assembler_Zoltan_test(argc, argv);
    return 0;
}































