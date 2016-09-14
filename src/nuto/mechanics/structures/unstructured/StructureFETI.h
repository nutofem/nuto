
#pragma once

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Sparse>

#include <cstring>
#include <fstream>


namespace NuTo
{


//! @author Philip Huschke
//! @date September 2016
//! @brief Class for FETI methods
class StructureFETI: public Structure
{
private:
    struct Node
    {
        Eigen::Vector3d mCoordinates;
        int mId;
    };

    struct Element
    {
        int mId;
        std::vector<int> mNodeIds;
    };

    struct Boundary
    {
        std::map<int, int> mNodeIdsMap;
        std::vector<int> mValues;
    };

    struct Interface
    {
        std::map<int, int> mNodeIdsMap;
        int mValue;
    };

    using NodeList      = std::vector<Node>;
    using ElementList   = std::vector<Element>;
    using BoundaryList  = std::vector<Boundary>;
    using InterfaceList = std::vector<Interface>;

    NodeList        ReadNodeData        (std::ifstream& file);
    ElementList     ReadElementData     (std::ifstream& file);
    BoundaryList    ReadBoundaryData    (std::ifstream& file);
    InterfaceList   ReadInterfaceData   (std::ifstream& file);


public:
    //! @brief Constructor
    //! @param rDimension   Structural dimension (1,2 or 3)
    //! @param rMeshFile    mesh file
    StructureFETI(int rDimension, std::string rMeshFile);

    const int mRank;
    const int mNumProcesses;

    void FindKeywordInFile(std::ifstream &file, std::string keyword);

    const Eigen::SparseMatrix<double>& GetConnectivityMatrix() {return mConnectivityMatrix;}
    Eigen::MatrixXd GetRigidBodyModes();
    Eigen::MatrixXd GetInterfaceRigidBodyModes();
    Eigen::SparseMatrix<double> &AssembleStiffnessMatrix();
    void AssembleRigidBodyModes();
    void AssembleConnectivityMatrix();
    Eigen::MatrixXd             mRigidBodyModes;
protected:
    void ImportMesh(std::string rFileName);


    Eigen::SparseMatrix<double> mConnectivityMatrix;

    Eigen::MatrixXd             mInterfaceRigidBodyModes;
    Eigen::MatrixXd             mProjectionMatrix;
    NodeList                    mNodes;
    ElementList                 mElements;
    BoundaryList                mBoundaries;
    InterfaceList               mInterfaces;

};
} //namespace NuTo

