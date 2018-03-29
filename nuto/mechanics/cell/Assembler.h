#pragma once
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/dofs/DofVector.h"
#include "nuto/mechanics/dofs/DofMatrix.h"
#include "nuto/mechanics/dofs/DofMatrixSparse.h"

namespace NuTo
{

//! Assembles an internal NuTo::DofVector<double> from arbitrary contributions and provides access to it
class VectorAssembler
{
public:
    //! ctor
    //! @param sizes number of global dofs for each dof type
    VectorAssembler(DofContainer<int> sizes = {});

    //! resizes the internal DofVector
    //! @param sizes number of global dofs for each dof type
    void Resize(DofContainer<int> sizes);

    //! adds an arbitrary contribution `v` to the internal DofVector
    //! @param v contribution, e.g. a local element/cell vector
    //! @param numbering mapping from local (0...numbering.size()) to global (0...sizes)
    //! @param dofTypes specific dof types to assemble - an empty value will assemble all available dof types from the
    //! numbering
    void Add(const DofVector<double>& v, const DofVector<int>& numbering, std::vector<DofType> dofTypes = {});

    //! sets the entries of the internal DofVector to zero
    void SetZero();

    //! getter
    //! @return current state of the internal DofVector
    const DofVector<double>& Get() const;

private:
    DofVector<double> mVector;
};


//! Assembles an internal NuTo::DofMatrixSparse<double> from arbitrary contributions and provides access to it
//! @remark an interal flag distinguishes from the first assembly (into triplet lists) and any other assemblies
//! (directly into the existsing nonzero entries)
class MatrixAssembler
{
public:
    //! ctor
    //! @param sizes number of global dofs for each dof type
    MatrixAssembler(DofContainer<int> sizes = {});

    //! resizes the internal DofMatrixSparse
    //! @param sizes number of global dofs for each dof type
    void Resize(DofContainer<int> sizes);

    //! adds an arbitrary contribution `m` to the internal DofMatrixSparse
    //! @param m contribution, e.g. a local element/cell matrix
    //! @param numbering mapping from local (0...numbering.size()) to global (0...sizes)
    //! @param dofTypes specific dof types to assemble - an empty value will assemble all available dof types from the
    //! numbering
    void Add(const DofMatrix<double>& m, const DofVector<int>& numbering, std::vector<DofType> dofTypes = {});

    //! adds an arbitrary contribution `m` to `rSparseMatrix`
    //! @param m contribution, e.g. a local element/cell matrix
    //! @param numbering mapping from local (0...numbering.size()) to global (0...sizes)
    //! @param dofTypes specific dof types to assemble - an empty value will assemble all available dof types from the
    //! numbering
    static void Add(DofMatrixSparse<double>& rSparseMatrix, const DofMatrix<double>& m, const DofVector<int>& numbering,
                    std::vector<DofType> dofTypes = {});

    //! sets the entries of the internal DofMatrixSparse to zero but keeps the nonzeros
    void SetZero();

    //! transforms the internal triplet lists into an DofMatrixSparse and drops the triplet lists
    //! @remark calling this method indicates that the position of the nonzeros are now fixed
    void Finish();

    //! @return current state of the internal DofMatrixSparse
    const DofMatrixSparse<double>& Get() const;

private:
    DofMatrixContainer<std::vector<Eigen::Triplet<double>>> mTriplets;
    DofMatrixSparse<double> mMatrix;
    bool mFinished = false;
};


} /* NuTo */
