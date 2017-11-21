#include <mpi.h>

#include <boost/mpi.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include "json.hpp"
#include "base/Exception.h"

#include "mechanics/feti/StructureFeti.h"
#include "mechanics/structures/unstructured/Structure.h"

#include <Epetra_MpiComm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Export.h>
#include <Epetra_LinearProblem.h>




class StructureMesh: public NuTo::StructureFeti
{
private:
    struct serializeData{
        std::vector<int> completeNodeIDs;
        std::vector<std::vector<double>> completeNodeCoords;
        std::vector<std::vector<int>> completeElementNodeIDs;
        std::vector<int> completeDisplacementIDs;
        std::vector<double> completeDisplacements;
        std::vector<std::vector<int>> completeNode2Dofs;

        template<typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::make_nvp;
            ar & make_nvp("completeNodeIDs", completeNodeIDs);
            ar & make_nvp("completeNodeCoords", completeNodeCoords);
            ar & make_nvp("completeElementNodeIDs", completeElementNodeIDs);
            ar & make_nvp("completeDisplacementIDs", completeDisplacementIDs);
            ar & make_nvp("completeDisplacements", completeDisplacements);
            ar & make_nvp("completeNode2Dofs", completeNode2Dofs);
        }
    };


    struct dofData{
        int numberDofs;
        int numberActiveDofs;
        int numberMasterDofs;
        int numberMasterActiveDofs;
        int numberMasterDependentDofs;

        template<typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::make_nvp;
            ar & make_nvp("numberDofs", numberDofs);
            ar & make_nvp("numberActiveDofs", numberActiveDofs);
            ar & make_nvp("numberMasterDofs", numberMasterDofs);
            ar & make_nvp("numberMasterActiveDofs", numberMasterActiveDofs);
            ar & make_nvp("numberMasterDependentDofs", numberMasterDependentDofs);
        }
    };




    std::vector<std::vector<double>> sortNodeCoords(std::vector<int> rNodeIDs, std::vector<std::vector<double>> rNodeCoords);


    std::vector<double> sortDisplacements(std::vector<int> rDisplacementIDs, std::vector<double> rDisplacements);


    serializeData serializeSolutionParticular(std::vector<double> rSolution, std::vector<int> rSolutionIDs, std::vector<std::vector<int>> rNode2Dofs, int rNumProc);


public:
    struct dofNode {
        Node nodeInfo;
        std::map<NuTo::Node::eDof, std::vector<int>> dofIDs;
        std::map<NuTo::Node::eDof, std::vector<int>> activeDofIDs;
        std::map<NuTo::Node::eDof, std::vector<int>> dependentDofIDs;
        int masterDomain;

        template<typename Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
            using boost::serialization::make_nvp;
            ar & make_nvp("nodeInfo", nodeInfo);
            ar & make_nvp("dofIDs", dofIDs);
            ar & make_nvp("activeDofIDs", activeDofIDs);
            ar & make_nvp("dependentDofIDs", dependentDofIDs);
            ar & make_nvp("masterDoamin", masterDomain);
        }
    };





    StructureMesh(int rDim) : NuTo::StructureFeti(rDim) {}

//    void importMyMeshJson(std::string rFileName, const int interpolationTypeId);
    void importMyMeshJson(std::string rFileName);


    std::vector<std::vector<int>> map2Vector(std::map<int, std::vector<int>> rMap);

    std::vector<std::vector<int>> map2Vector(std::map<int, Eigen::VectorXi> rMap);


    std::vector<int> map2Vector(std::map<int, int> rMap);

    std::vector<int> map2ValueVector(std::map<int, int> rMap);

    void visualizeSerializedParticularSolution(std::vector<double> rSolution, std::vector<int> rSolutionIDs, std::vector<std::vector<int>> rNode2Dof, std::string rFileName, int rNumProc);


    void visualizeSolution(std::vector<double> rSolution, std::string rFileName);



    void generateNodeToDofMapping();


    void generateDofClassification();


    void gatherNodeToDofMapping_allProcesses(int rNumProc);


    void gatherNodeToMasterDofMapping(int rNumProc, int rRank);


    void gatherNodeToSlaveDofMapping(int rNumProc, int rRank);


    void gatherNodeToDofMapping(int rNumProc, int rRank);


    Epetra_MultiVector solveSystem(Epetra_CrsMatrix rA, Epetra_MultiVector rRhs, bool iterative = true, bool useAztecOO = true);


    Epetra_MultiVector solveSystem(Epetra_CrsMatrix rA, Epetra_MultiVector rLhs, Epetra_MultiVector rRhs, bool iterative = true, bool useAztecOO = true);


    std::map<int, Eigen::VectorXi> getNodeActiveDOFs(NuTo::Node::eDof rDofType, int rRank);

    std::map<int, Eigen::VectorXi> getNodeDOFs(NuTo::Node::eDof rDofType, int rRank);

    std::map<int, Eigen::VectorXi> getMyNodeDOFs(NuTo::Node::eDof rDofType);

    std::vector<int> getAllNodeIDs();

    std::vector<std::map<NuTo::Node::eDof, std::map<int,int>>> getLocalToGlobalDofMapping(){ return mLocalToGlobalDofMapping;}

    std::vector<std::map<NuTo::Node::eDof, std::map<int,int>>> getLocalToGlobalActiveDofMapping(){ return mLocalToGlobalActiveDofMapping;}

    std::vector<std::map<NuTo::Node::eDof, std::vector<int>>> getMasterGlobalDofIDs(){ return mMasterGlobalDofs;}

    std::vector<std::map<NuTo::Node::eDof, std::vector<int>>> getMasterGlobalActiveDofIDs(){ return mMasterGlobalActiveDofs;}

    std::vector<std::map<NuTo::Node::eDof, std::vector<int>>> getMasterGlobalDependentDofIDs(){ return mMasterGlobalDependentDofs;}

    std::vector<dofNode> getAllNodes(){ return mMyNodes;}

    std::map<NuTo::Node::eDof, std::map<int, int>> getMyLocalToGlobalDofMapping(){ return mMyLocalToGlobalDofMapping;}

    std::map<NuTo::Node::eDof, std::map<int, int>> getMyLocalToGlobalActiveDofMapping(){ return mMyLocalToGlobalActiveDofMapping;}

    std::map<NuTo::Node::eDof, std::vector<int>> getMyMasterGlobalDofIDs(){ return mMyMasterGlobalDofs;}

    std::map<NuTo::Node::eDof, std::vector<int>> getMyMasterGlobalActiveDofIDs(){ return mMyMasterGlobalActiveDofs;}

    std::map<NuTo::Node::eDof, std::vector<int>> getMyMasterGlobalDependentDofIDs(){ return mMyMasterGlobalDependentDofs;}


protected:
    std::vector<dofNode> mMyNodes;
    std::vector<dofNode> mMyMasterNodes;
    std::vector<dofNode> mMySlaveNodes;

    std::vector<std::map<NuTo::Node::eDof, std::map<int, int>>> mLocalToGlobalDofMapping;
    std::vector<std::map<NuTo::Node::eDof, std::map<int, int>>> mLocalToGlobalActiveDofMapping;
    std::vector<std::map<NuTo::Node::eDof, std::vector<int>>> mMasterGlobalDofs;
    std::vector<std::map<NuTo::Node::eDof, std::vector<int>>> mMasterGlobalActiveDofs;
    std::vector<std::map<NuTo::Node::eDof, std::vector<int>>> mMasterGlobalDependentDofs;

    std::map<NuTo::Node::eDof, dofData> mMyDofData;
    std::map<NuTo::Node::eDof, std::map<int, int>> mMyLocalToGlobalDofMapping;
    std::map<NuTo::Node::eDof, std::map<int, int>> mMyLocalToGlobalActiveDofMapping;
    std::map<NuTo::Node::eDof, std::vector<int>> mMyMasterGlobalDofs;
    std::map<NuTo::Node::eDof, std::vector<int>> mMyMasterGlobalActiveDofs;
    std::map<NuTo::Node::eDof, std::vector<int>> mMyMasterGlobalDependentDofs;

};

