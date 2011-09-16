// $Id$

#ifndef STRUCTUREMULTISCALE_H
#define STRUCTUREMULTISCALE_H

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/export.hpp>
#include <boost/serialization/array.hpp>
#else
#include <boost/array.hpp>
#endif //Serialize
#include "boost/filesystem.hpp"

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveEngineeringStressStrain.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress2D.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"


namespace NuTo
{
template <class T> class FullMatrix;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard class for irregular (unstructured) structures, which are used as fine scale models at an integration point
class StructureMultiscale : public Structure
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
    //just for test purpose
    friend class ConstraintNonlinearGlobalCrackAngle2D;
public:
    using StructureBase::NewtonRaphson;
    //! @brief constructor
    //! @param mDimension  Structural dimension (1,2 or 3)
    StructureMultiscale(int mDimension);

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif //SWIG

    //! @brief ... save the object to a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Save (const std::string &filename, std::string rType )const;

    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    void Restore (const std::string &filename, std::string rType );
#endif // ENABLE_SERIALIZATION

    //! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const
    {
        return std::string("StructureMultiscale");
    }

    //************ Info routine         ***************
    //**  defined in structures/StructureMultiscale.cpp *********
    //*************************************************
    //! @brief ... Info routine that prints general information about the object (detail according to verbose level)
    void Info()const;

#ifndef SWIG
    //! @brief ... calculate the displacement based on the homogeneous strain
    //! @param rCoordinates ... coordinates of the point
    //! @param rCoarseDisplacements ... return array of displacements
    void GetDisplacementsEpsilonHom2D(double rCoordinates[2], double rDisplacements[2], const boost::array<double,2>& rCenter)const;

    //! @brief ... calculate the displacement based on the crack opening
    //! @param rCoordinates ... coordinates of the point
    //! @param rCoarseDisplacements ... return array of displacements
    void GetDisplacementsCrack2D(double rCoordinates[2], double rDisplacements[2])const;

    //! @brief derivative of displacement with respect to homogeneous strain
    //! @param rdX_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
    //! @param rdY_dEpsilonHom[3] return value, derivative of x-displacement with respect to homogeneous strain (exx, eyy, gxy)
    void GetdDisplacementdEpsilonHom(double rCoordinates[2], double rdX_dEpsilonHom[3], double rdY_dEpsilonHom[3], const boost::array<double,2>& rCenter)const;

    //! @brief derivative of displacement with respect to discontinuity (crack opening)
    //! @param rdX_dCrackOpening[2] return value, derivative of x-displacement with respect to crack opening (ux, uy)
    //! @param rdY_dCrackOpening[2] return value, derivative of x-displacement with respect to crack opening (ux, uy)
    void GetdDisplacementdCrackOpening(double rCoordinates[2], double rdX_dCrackOpening[2], double rdY_dCrackOpening[2])const;

#endif //SWIG

    //! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
    //! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
    void TransformMultiscaleNodes(int rGroupBoundaryNodes, int rGroupNodes, int rGroupElements, boost::array<double,2>& rCenter, double& rArea, bool rCrackedDomain);

    //! @brief ... the boundary nodes were transformed from pure displacement type nodes to multiscale nodes
    //! the displacements are decomposed into a local displacement field and a global homogeneous/crack displacement
    void TransformMultiscaleNodes();

#ifndef SWIG
    //! @brief numbers non standard DOFs' e.g. in StructureMultiscale, for standard structures this routine is empty
    void NumberAdditionalGlobalDofs();

    //! @brief ...merge additional dof values
    void NodeMergeAdditionalGlobalDofValues(const NuTo::FullMatrix<double>& rActiveDofValues, const NuTo::FullMatrix<double>& rDependentDofValues);

    //! @brief extract dof values additional dof values
    //! @param rActiveDofValues ... vector of global active dof values (ordering according to global dofs, size is number of active dofs)
    //! @param rDependentDofValues ... vector of global dependent dof values (ordering according to (global dofs) - (number of active dofs), size is (total number of dofs) - (number of active dofs))
    void NodeExtractAdditionalGlobalDofValues(NuTo::FullMatrix<double>& rActiveDofValues, NuTo::FullMatrix<double>& rDependentDofValues) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatrices0General(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKJ, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK) const;

    //! @brief ... based on the global dofs build submatrices of the global coefficent matrix0
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    NuTo::Error::eError BuildGlobalCoefficientSubMatrices0Symmetric(NuTo::SparseMatrix<double>& rMatrixJJ, NuTo::SparseMatrix<double>& rMatrixJK, NuTo::SparseMatrix<double>& rMatrixKK) const;

    //! @brief ... for a single element matrix, assemble the global matrix, consider here additionally the global degrees of freedom from StructureMultiscale (alpha, ux,uy) and eventually the macroscopic strain
    //! @param elementMatrix ...stiffness matrix of the element
    //! @param elementMatrixGlobalDofsRow ...corresponding row DOFs
    //! @param elementMatrixGlobalDofsColumn ...corresponding column DOFs
    //! @param elementVector ...gradient matrix of the element (internal force)
    //! @param mappingDofMultiscaleNode ...mapping of all Dofs to Multiscale dofs (either -1 for standard dofs or the entry in dDOF and dDOF2, where additional infos are stored)
    //! @param dDOF ... derivatives of multiscale dofs with respect to alpha, ux, uy, ehomxx, ehomyy, gammahomxy [0..6] scaled derivatives [7..12]
    //! @param dDOF2 ...for each dof, the corresponding second order derivative (alpha^2, alpha ux, alpha uy)
    //! @param rMatrixJJ ... submatrix jj (number of active dof x number of active dof)
    //! @param rMatrixJK ... submatrix jk (number of active dof x number of dependent dof)
    //! @param rMatrixKJ ... submatrix kj (number of dependent dof x number of active dof)
    //! @param rMatrixKK ... submatrix kk (number of dependent dof x number of dependent dof)
    //! @param rCalcMatrixKJ_KK ...true, if rMatrixKJ and rMatrixKK should be calculated (rMatrixJJ and rMatrixJK) are always calculated
    //! @param rCalcMatrixM ...true, if all matrices corresponding to M should be calculated
    void AddElementMatrixToGlobalSubMatricesGeneral(
            NuTo::FullMatrix<double>& elementMatrix,
            std::vector<int>& elementMatrixGlobalDofsRow,
            std::vector<int>& elementMatrixGlobalDofsColumn,
            std::vector<int>& mappingDofMultiscaleNode,
            std::vector<boost::array<double,5> >& rDOF,
            NuTo::SparseMatrix<double>* rMatrixJJ,
            NuTo::SparseMatrix<double>* rMatrixJK,
            NuTo::SparseMatrix<double>* rMatrixKJ,
            NuTo::SparseMatrix<double>* rMatrixKK,
            bool rCalcMatrixKJ_KK)const;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    NuTo::Error::eError BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector) const;

    //! @brief ... based on the global dofs build sub-vectors of the global internal potential gradient
    //! @param rActiveDofGradientVector ... global internal potential gradient which corresponds to the active dofs
    //! @param rDependentDofGradientVector ... global internal potential gradient which corresponds to the dependent dofs
    //! @param rEpsilonMGradientVector ... global internal potential gradient which corresponds to the global macroscopic strain (derivative is equal to homogeneous strain)
    NuTo::Error::eError BuildGlobalGradientInternalPotentialSubVectors(NuTo::FullMatrix<double>& rActiveDofGradientVector, NuTo::FullMatrix<double>& rDependentDofGradientVector, NuTo::FullMatrix<double>* rEpsilonMGradientVector) const;

    void CalculateCohesiveForce(NuTo::FullMatrix<double>& rCohesiveForce) const;

    //! @brief calculate the distance of a point to the crack
    //! @param rCoordinates[2] coordinates of the point
    //! @return distance to crack
    double CalculateDistanceToCrack2D(double rCoordinates[2])const;

    //! @brief calculate the derivative of the displacements at the nodes with respect to homogeneous strain, crack opening and crack orientation
    //! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
    //! @param rDOF return value, for each dof, the corresponding derivatives ehomxx, ehomyy, gammahomxy, ux, uy,
    void CalculatedDispdGlobalDofs(std::vector<int>& rMappingDofMultiscaleNode, std::vector<boost::array<double,5> >& rDOF)const;

    //! @brief calculate the derivative of the displacements at the nodes with respect to crack opening without considering this to be balanced by epsilon_hom
    //! @param rMappingDofMultiscaleNode return value, for each dof, the corresponding entry in the rDOF vector, for nonmultiscale dofs, there is a -1
    //! @param rDOF return value, for each dof, the corresponding derivatives (ux, uy)
    void CalculatedDispdCrackOpening(std::vector<int>& rMappingDofMultiscaleNode,std::vector<boost::array<double,2> >& rDOF)const;

    //! @brief set the total current strain
    void SetTotalEngineeringStrain(const EngineeringStrain2D& rTotalEngineeringStrainConstraint);

    //! @brief set the total current strain as a constraint, almost identical to totalStrain
    void SetTotalEngineeringStrainConstraint(const EngineeringStrain2D& rTotalEngineeringStrainConstraint);

    //! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
    void SetDeltaTotalEngineeringStrain(const EngineeringStrain2D& rDeltaTotalEngineeringStrain);

    //! @brief set the maximum total strain (used in the automatic increment of the Newton Raphson iteration multiplied by the load factor to give the totalEngineeringStrain)
    void SetPrevTotalEngineeringStrain(const EngineeringStrain2D& rPrevTotalEngineeringStrain);

    //! @brief returns the total strain
    const NuTo::EngineeringStrain2D& GetTotalEngineeringStrain()const;

    //! @brief returns the total strain
    const NuTo::EngineeringStrain2D& GetTotalEngineeringStrainConstraint()const;

    //! @brief just for testing
    void SetGlobalCrackOpening(boost::array<double,2> ptr)const
    {
        const_cast<StructureMultiscale*>(this)->mCrackOpening[0] = ptr[0];
        const_cast<StructureMultiscale*>(this)->mCrackOpening[1] = ptr[1];
    }

    //! @brief set crack angle
    void SetCrackAngle(double alpha)
    {
        mCrackAngle = alpha;
    }

    //! @brief return the total strain
    const NuTo::EngineeringStrain2D& GetHomogeneousEngineeringStrain()const;

    //! @brief set the homgeneous strain (just for test purpose)
    void SetHomogeneousEngineeringStrain(NuTo::EngineeringStrain2D rStrain);

    //! @brief renumbers the global dofs in the structure after
    void ReNumberAdditionalGlobalDofs(std::vector<int>& rMappingInitialToNewOrdering);
#endif
    //! @brief returns the constraint equation for the total strain
    int GetConstraintTotalStrain()const
    {
        return mConstraintTotalStrain;
    }

    //! @brief returns the constraint equation for the tangential crack opening
    int GetConstraintCrackOpeningTangential()const
    {
        return mConstraintTangentialCrackOpening;
    }

    //! @brief returns the constraint equation for the normal crack opening
    int GetConstraintCrackOpeningNormal()const
    {
        return mConstraintNormalCrackOpening;
    }

    void SetlCoarseScaleCrack(double rlCoarseScaleCrack)
    {
        mlCoarseScaleCrack = rlCoarseScaleCrack;
    }

    double GetlCoarseScaleCrack()const
    {
        return mlCoarseScaleCrack;
    }

    double GetAreaFineScale()const
    {
        return mFineScaleArea;
    }

    double GetCoarseScaleArea()const
    {
        return mCoarseScaleArea;
    }

    void SetCoarseScaleArea(double rCoarseScaleArea)
    {
        mCoarseScaleArea=rCoarseScaleArea;
    }

    double GetCrackAngle()const
    {
        return mCrackAngle;
    }

    const boost::array<int,2>& GetDofGlobalCrackOpening2D()const
    {
        return mDOFCrackOpening;
    }

    const boost::array<double,2>& GetGlobalCrackOpening2D()const
    {
        return mCrackOpening;
    }
    const boost::array<int,3>& GetDofGlobalTotalStrain2D()const
    {
        return mDOFGlobalTotalStrain;
    }

    const boost::array<double,3>& GetPeriodicBoundaryDisplacements()const
    {
    	return mPeriodicBoundaryDisplacements;
    }

    const boost::array<int,3>& GetDOFPeriodicBoundaryDisplacements()const
    {
    	return mDOFPeriodicBoundaryDisplacements;
    }

    void SetGroupBoundaryNodesElements(int rGroupIdBoundaryNodes, int rGroupIdNodes, int rGroupIdElements)
    {
        mGroupBoundaryNodesHomogeneous = rGroupIdBoundaryNodes;
        mGroupNodesHomogeneous = rGroupIdNodes;
        mGroupElementsHomogeneous = rGroupIdElements;
        mBoundaryNodesElementsAssigned = true;
    }

    int GetGroupElementsDamage()const
    {
    	return mGroupElementsDamage;
    }

    int GetGroupElementsHomogeneous()const
    {
    	return mGroupElementsHomogeneous;
    }

    int GetGroupBoundaryNodesDamage()const
    {
    	return mGroupBoundaryNodesDamage;
    }

    int GetGroupBoundaryNodesHomogeneous()const
    {
    	return mGroupBoundaryNodesHomogeneous;
    }

    int GetGroupNodesDamage()const
    {
    	return mGroupNodesDamage;
    }

    int GetGroupNodesHomogeneous()const
    {
    	return mGroupNodesHomogeneous;
    }

    double GetScalingFactorDamage()const
    {
    	assert(mlCoarseScaleCrack/mlFineScaleCrack>0);
    	assert(mGroupElementsDamage!=-1);
    	//std::cout << "scaling Factor Dam " << mlCoarseScaleCrack/mlFineScaleCrack << "\n";
    	//std::cout << "mlFineScaleCrack " <<mlFineScaleCrack << "\n";
    	//std::cout << "mlCoarseScaleCrack " << mlCoarseScaleCrack << "\n";
        return mlCoarseScaleCrack/mlFineScaleCrack;
    }

    double GetScalingFactorHomogeneous()const
    {
/*    	std::cout << "mCoarseScaleArea " << mCoarseScaleArea << "\n";
    	std::cout << "mlCoarseScale " << mlCoarseScale << "\n";
    	std::cout << "mFineScaleAreaDamage " << mFineScaleAreaDamage << "\n";
    	std::cout << "mlFineScaleDamage " << mlFineScaleDamage << "\n";
    	std::cout << "mFineScaleAreaHomogeneous " << mFineScaleAreaHomogeneous << "\n";
*/
        if (mGroupElementsDamage==-1)
        {
        	//std::cout << "scaling Factor Hom " << mCoarseScaleArea/mFineScaleArea << "\n";
        	assert(mCoarseScaleArea/mFineScaleArea>0);
        	return mCoarseScaleArea/mFineScaleArea;
        }
        else
        {
        	//std::cout << "scaling Factor Hom " << (mCoarseScaleArea-mlCoarseScaleCrack*mFineScaleArea/mlFineScaleCrack)/(mFineScaleArea) << "\n";
        	assert((mCoarseScaleArea-mlCoarseScaleCrack*mFineScaleArea/mlFineScaleCrack)/(mFineScaleArea)>0);
    	    return (mCoarseScaleArea-mlCoarseScaleCrack*mFineScaleArea/mlFineScaleCrack)/(mFineScaleArea);
        }
    }

    const std::string& GetIPName()const
    {
        return mIPName;
    }

    void SetIPName(std::string rIPName)
    {
        mIPName = rIPName;
    }

    void SetCenterMacro(double rCenterMacro[2])
    {
    	mCenterMacro[0] = rCenterMacro[0];
    	mCenterMacro[1] = rCenterMacro[1];
    }

    void SetShiftCenterDamage(double rShiftCenterDamage[2])
    {
    	mShiftCenterDamage[0] = rShiftCenterDamage[0];
    	mShiftCenterDamage[1] = rShiftCenterDamage[1];
    }

    void SetCrackOpening(NuTo::FullMatrix<double>& crackOpening);

    void GetShiftedCenterDamage(boost::array<double,2> rShiftedCenterDamage)const
    {
    	rShiftedCenterDamage[0] = mCenterDamage[0] + mShiftCenterDamage[0];
    	rShiftedCenterDamage[1] = mCenterDamage[1] + mShiftCenterDamage[1];
    }

    //calculate from the existing crack opening and orientation the cracking strain and the homogeneous strain
    //for test purpose in public section
    void CalculateHomogeneousEngineeringStrain();

    void SetCrackTransitionRadius(double rCrackTransitionRadius)
    {
        mCrackTransitionRadius = rCrackTransitionRadius;
    }

    //! @brief add a constraint equation for the crack opening (normal crack opening non negativ)
    //! @parameter rPenaltyStiffness penalty stiffness for augmented Lagrangian
    //! @return id of the constraint
    int CreateConstraintLagrangeGlobalCrackOpeningNormal(double rPenaltyStiffness);

    //! @brief add a constraint equation for the total strain
    //! @parameter rStrain applied strain (rhs)
    //! @return id of the constraint
    int CreateConstraintLinearGlobalTotalStrain();

    //! @brief add a linear constraint equation for the crack opening
    //! @return id of the constraint
    int CreateConstraintLinearGlobalCrackOpening(double rRHS, const NuTo::FullMatrix<double>& rDirection);

    //! @brief add a linear constraint equation for the crack opening in tangential direction
    //! @return id of the constraint
    int CreateConstraintLinearGlobalCrackOpeningTangential(double rRHS);

    //! @brief add a linear constraint equation for the crack opening in tangential direction
    //! @return id of the constraint
    int CreateConstraintLinearGlobalCrackOpeningNormal(double rRHS);

    //! @brief add a linear constraint equation for the additional shape functions depscribing the fluctuation boundary displacements
    //! @return id of the constraint
    int CreateConstraintLinearPeriodicBoundaryShapeFunctions(int rShapeFunction, double rRHS);

    //! @brief set constraint for fine scale fluctuations on the boundary as a linear combination of the periodic bc with exx, eyy, gxy
    //! @parameter rHomogeneousDomain true for the homogeneous domain, false for the damage domain
    void CreateConstraintLinearFineScaleDisplacementsUsingAddShapeFunctions(bool rHomogeneousDomain);

    //! @brief add constraints for the finescale displacement of the boundary nodes
    //! @parameter rDispX .. displacement in x-direction
    //! @parameter rDispY .. displacement in y-direction
    //! @parameter rHomogeneousDomain .. true for the homogeneous domain, false for the damage domain
    void CreateConstraintLinearFineScaleDisplacements(double rDispX, double rDispY, bool rHomogeneousDomain);

    //! @brief this routine is only relevant for the multiscale model, since an update on the fine scale should only be performed
    //for an update on the coarse scale
    //as a consequence, in an iterative solution with updates in between the initial state has to be restored after leaving the routine
    //this routine saves the current state before an update in the Newton Raphson iteration is performed
    //this only happens for more than one load step (either prescibed or with automatic load control)
    //if not already saved, the structure is saved to an internal stream
    void SaveStructure(std::stringstream& rSaveStringStream)const;

    //restores the structure previously saved into the string
    void RestoreStructure(std::stringstream& rSaveStringStream);

    //! @brief set the load factor (load or displacement control) overload this function to use Newton Raphson
    //! @param load factor
    void SetLoadFactor(double rLoadFactor);

    //! @brief do a postprocessing step after each converged load step (for Newton Raphson iteration) overload this function to use Newton Raphson
    void PostProcessDataAfterConvergence(int rLoadStep, int rNumNewtonIterations, double rLoadFactor, double rDeltaLoadFactor, double rResidual)const;

    //! @brief do a postprocessing step after each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
    void PostProcessDataAfterLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, const NuTo::FullMatrix<double>& rResidualVector)const;

    //! @brief do a postprocessing step in each line search within the load step(for Newton Raphson iteration) overload this function to use Newton Raphson
    void PostProcessDataInLineSearch(int rLoadStep, int rNewtonIteration, double rLineSearchFactor, double rLoadFactor, double rResidual, double rPrevResidual)const;

    //! @brief initialize some stuff before a new load step (e.g. output directories for visualization, if required)
    void InitBeforeNewLoadStep(int rLoadStep);

    double GetlFineScaleCrack()const
    {
        return mlFineScaleCrack;
    }

    void SetlFineScaleCrack(double rlFineScaleCrack)
    {
        mlFineScaleCrack = rlFineScaleCrack;
    }

    double GetlFineScale()const
    {
        return mlFineScale;
    }

    const boost::array<double,2>& GetCenterDamage()const
    {
        return mCenterDamage;
    }

    const boost::array<double,2>& GetCenterHomogeneous()const
    {
        return mCenterHomogeneous;
    }

    //! @brief calculates the total energy of the system
    //! @return total energy
    double ElementTotalGetTotalEnergy()const;

    //! @brief calculates the total energy of the system
    //! @return total energy
    double ElementGroupGetTotalEnergy(int rGroupId)const;

    //! @brief sets the result directory where the results are written to
    void SetResultDirectory(std::string rResultDirectory);

    //! @brief gets the result directory where the results are written to
    std::string GetResultDirectory()const
    {
    	return mResultDirectory;
    }

    //! @brief sets the result file for the converged solution where the results are written to
    void SetResultLoadStepMacro(std::string rLoadStepMacro);

    //! @brief gets the result file for the converged solution where the results are written to
    std::string GetResultLoadStepMacro()const
    {
    	return mLoadStepMacro;
    }
    //! @brief returns the file where the log is written to for each integration point
    std::string GetFileForLog()const
    {
    	boost::filesystem::path resultDir(mResultDirectory);
    	resultDir /= mIPName;
    	resultDir /= std::string("output_") + mLoadStepMacro + std::string(".log");
    	return resultDir.string();
    }

    //! @brief returns the file where the visualization is written to for each integration point
    std::string GetFileForVisualizationAfterLineSearch(int rLoadStep, int rNewtonIteration)const
    {
    	boost::filesystem::path resultDir(mResultDirectory);
    	resultDir /= mIPName;
    	std::stringstream ssLoadStep;
    	ssLoadStep << rLoadStep;
    	std::stringstream ssNewtonIteration;
    	ssNewtonIteration << rNewtonIteration;
    	resultDir /= std::string("Structure_") + mLoadStepMacro + "_" + ssLoadStep.str() +"_" + ssNewtonIteration.str() + std::string(".vtk");
    	return resultDir.string();
    }

    //! @brief is only true for structure used as multiscale (structure in a structure)
    virtual bool IsMultiscaleStructure()const
    {
    	return true;
    }

    //! @brief is only true for structure used as multiscale (structure in a structure)
    void ScaleCoordinates(double rCoordinates[3])const;

    //! @brief sets the center for scaling the coordinates either to the center of the damage model(true) or the homogeneous model (false)
    void SetCenterScalingToDamage(bool rScaleWRTDamageCenter)
    {
    	mScaleWRTDamageCenter = rScaleWRTDamageCenter;
    }

#ifdef ENABLE_VISUALIZE
    // visualizes the crack of the fine scale model
    void VisualizeCrack(VisualizeUnstructuredGrid& rVisualize)const;
#endif //ENABLE_VISUALIZE

/*    //! @brief performs a Newton Raphson iteration (displacement and/or load control)
    //! @parameters rSaveStructureBeforeUpdate if set to true, save the structure (done in a separate routine to be implemented by the user) before an update is performed
    //!             be careful, store it only once, although the routine is called before every update
    void NewtonRaphson(bool rSaveStructureBeforeUpdate,
            std::stringstream& rSaveStringStream,
            bool& rIsSaved);
*/
    void CalculateStiffness(NuTo::FullMatrix<double>& rStiffness, bool rPeriodic);

    void CalculatePeriodicBoundaryShapeFunctions(double rDeltaStrain);

    inline void SetScalingFactorCrackOpening(double rScalingFactorCrackOpening)
    {
    	mScalingFactorCrackOpening = rScalingFactorCrackOpening;
    }

    inline double GetScalingFactorCrackOpening()const
    {
    	return mScalingFactorCrackOpening;
    }

    inline void SetScalingFactorEpsilon(double rScalingFactorEpsilon)
    {
    	mScalingFactorEpsilon = rScalingFactorEpsilon;
    }

    inline double GetScalingFactorEpsilon()const
    {
    	return mScalingFactorEpsilon;
    }

    //! @brief delete the constraint for the tangential crack opening
    void ConstraintDeleteTangentialCrackOpening();

    //! @brief delete the constraint for the normal crack opening
	void ConstraintDeleteNormalCrackOpening();

#ifndef SWIG
	//! @brief set periodic boundary conditions for a 2D structure
	//! @parameter rGroupBoundaryNodes ... boundary nodes
	//! @parameter rStrain ... strain
	int ConstraintLinearSetFineScaleDisplacementsPeriodicNodeGroup(int rGroupBoundaryNodes, const EngineeringStrain2D& rStrain);
#endif

	//! @brief performs a Newton Raphson iteration (displacement and/or load control)
    //! @parameters rSaveStructureBeforeUpdate if set to true, save the structure (done in a separate routine to be implemented by the user) before an update is performed
    //!             be careful, store it only once
//	NuTo::StructureEnum::eErrorCode NewtonRaphson(bool rSaveStructureBeforeUpdate,
//            std::stringstream& rSaveStringStream,
//            bool& rIsSaved, bool rInitialStateInEquilibrium);

    //! @brief the global matrix is checked, not the element matrices
    //! @return false, if stiffness is not correct
    bool CheckStiffness();

    //! @brief calculates the crack angle for elastic solutions (initial value, no scaling with previous crack angle)
    //! @return alpha crack angle in the range 0..Pi (return value)
    double CalculateCrackAnglePrincipalStrain(const EngineeringStrain2D& rStrain)const;

	//! @brief calculates the length of a crack in the fine scale
    //! be careful, that the length, square/round, the angle and the offset of the finescale have to be set before
    void CalculateAndSetCrackLengthFineScale();

    //! @brief translates the homogeneous domain to the right and recreates the groups for boundary nodes, nodes and elements for the damaged domain
    void CreateDamageDomainFromHomogeneousDomain();

protected:
    //! @brief ... standard constructor just for the serialization routine
    StructureMultiscale();

    //! @brief Calculate the derivate of the homogeneous strain with respect to changes of the crack opening
    //! this is due to the constraint equation relating total strain, homogeneous strain and cracking strain
    //! @paramter rbHomU[0-2] for wrt ux [3-5] for wrt uy
    void GetdEpsilonHomdCrack(double rbHomU[6])const;

    double mCrackAngle;
    boost::array<double,2>  mShiftCenterDamage; //relative to center
    boost::array<double,2> mCrackOpening; //UT UN
    boost::array<int,2> mDOFCrackOpening;
    //! @brief this is the current total strain
    EngineeringStrain2D mEpsilonTot;
    //! @brief this is the strain related to the constrained, in general this is identical to mEpsilonTot
    //! but only before the load application, the stiffness is to be calculated with the old strain, but the
    //! constraint matrix with the new one
    EngineeringStrain2D mEpsilonTotConstraint;
    //! @brief dofs for the total strain (used in the stiffness calculation, otherwise this is always constrained)
    boost::array<int,3> mDOFGlobalTotalStrain;
    //! @brief dofs displacements of the boundary based on the solution with exx, eyy, and gxy
    boost::array<int,3> mDOFPeriodicBoundaryDisplacements;
    //! @brief corresponding displacements d=N(periodic shape functions of multiscale node)*mPeriodicBoundaryDisplacements
    boost::array<double,3> mPeriodicBoundaryDisplacements;
    //! @brief this is the current homogeneous part of the strain
    EngineeringStrain2D mEpsilonHom;
    //scaling factor for un/ut=mScalingFactorCrackOpening*DOFun/ut
    double mScalingFactorCrackOpening;
    //scaling factor for epsilonTot=mScalingFactorEpsilon*DOFepsilonTot
    double mScalingFactorEpsilon;
    //! @brief for length of the crack in the macroscopic element
    double mlCoarseScaleCrack;
    //! @brief parameter length of the fine scale crack in the damage model
    double mlFineScaleCrack;
    //! @brief parameter length of the fine scale (round=diameter, square=edgelength)
    double mlFineScale;
   //! @brief length to regularize the Heaviside function
    double mCrackTransitionRadius;
    //! @brief center of the model
    boost::array<double,2> mCenterDamage;
    //! @brief center of the model
    boost::array<double,2> mCenterHomogeneous;
    //! @brief center of the model (ip coordinates)
    boost::array<double,2> mCenterMacro;
    //! @brief translate and scale the coordinates for ip visualization either with respect to the damage center or the homogeneous center
    bool mScaleWRTDamageCenter;
    //! @brief area of the fine scale model
    // this has to be stored, since for a round fine scale model (especially coarse scale) the real area is substantially
    // lower than the area of the circle
    double mFineScaleArea;
    double mCoarseScaleArea;
    bool mSquareFineScaleModel;
    //! @brief thickness of the structure
    double mThickness;
    int mConstraintFineScaleDamageX,
        mConstraintFineScaleDamageY,
        mConstraintFineScaleHomogeneousX,
        mConstraintFineScaleHomogeneousY,
        mConstraintFineScalePeriodicDamage,
        mConstraintFineScalePeriodicHomogeneous,
        mConstraintNormalCrackOpening,
        mConstraintTangentialCrackOpening,
        mConstraintTotalStrain,
        mConstraintPeriodicBoundaryShapeFunctions[3];
    bool mBoundaryNodesElementsAssigned;
    int mGroupBoundaryNodesDamage;
    int mGroupBoundaryNodesHomogeneous;
    int mGroupNodesDamage;
    int mGroupNodesHomogeneous;
    int mGroupElementsDamage;
    int mGroupElementsHomogeneous;
    std::string mIPName;
    //result file for ips is stored in mResultDirectory/mIPName/mLoadStepMacro
    std::string mResultDirectory;
    //results for the whole structure are stored in mResultDirectory/mLoadStepMacro
    std::string mLoadStepMacro;
    // the following parameters are used for the Newton-Raphson iteration
    //! @brief this is the previous strain (used in the load application of the Newton procedure)
    EngineeringStrain2D mPrevEpsilonTot;
    //! @brief this is the delta strain (used in the load application of the Newton procedure)
    EngineeringStrain2D mDeltaEpsilonTot;
/*    //! @brief set to true, if during an update in NR the structure had to be saved, otherwise unchanged
    mutable bool mSavedToStringStream;
    //! @brief structure is saved to that string stream
    mutable std::stringstream mSaveStringStream;
*/
};
}
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::StructureMultiscale)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
#endif // StructureMultiscale_H
