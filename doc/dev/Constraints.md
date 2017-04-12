@page Constraints Handling of constraints in NuTo

# Theory

From the principle of virtual displacements, it follows for a body in equilibrium:
$$
\delta W_{int}^{(t+1)}= \delta W_{ext}^{(t+1)},
$$
where the variation of the internal energy $\delta W_{int}^{(t+1)}$ and the variation of the external energy $\delta W_{ext}^{(t+1)}$ is given by
\f{eqnarray*}{
\delta W_{int}^{(t+1)}&=\int_V\boldsymbol{\sigma}^T(\boldsymbol{d}^{(t+1)})\delta\boldsymbol{\varepsilon}^{(t+1)}\;dV\\
&=\int_V\left(\boldsymbol{\sigma}(\boldsymbol{d}^{(t)})+\Delta \boldsymbol\sigma\left(\boldsymbol{d}^{(t)}+\Delta \boldsymbol{d}^{(t+1)}\right)\right)^T\delta\boldsymbol\varepsilon^{(t+1)}\;dV\\
\delta W_{ext}^{(t+1)}&=\left(\boldsymbol F_{ext}^{(t)}+\Delta \boldsymbol{F}_{ext}^{(t+1)}\right)^T \delta \boldsymbol{d}^{(t+1)}.
\f}

Linearization of the system of equations gives
\f{align*}{
\int_V\left(\dfrac{\partial \sigma(\boldsymbol{d}^{(t)})}{\partial \varepsilon}\Delta\boldsymbol\varepsilon^{(t+1)}\right)^T\delta\boldsymbol\varepsilon^{(t+1)}\;dV=
\left(\Delta \boldsymbol F_{ext}^{(t+1)}\right)^T\delta \boldsymbol{d}^{(t+1)}-\int_V\left(\boldsymbol\sigma^{(t)}\right)^T\delta\boldsymbol\varepsilon^{(t+1)}\;dV+\left(\boldsymbol{F}_{ext}^{(t)}\right)^T\delta\boldsymbol d^{(t+1)},
\f}
where the last two summands on the right hand side vanish, if the previous step is in equilibrium. 
In the general formulation, the vector \f$\boldsymbol{d}\f$ includes all the nodal displacements. However, if kinematic boundary conditions have to be applied, the variations of the displacement state have to be kinematically compatible. In general, the linear constraint equations can be written as
\f{align*}{ 
\boldsymbol{T}_1^{(t+1)} \boldsymbol{d}_1^{(t+1)}+\boldsymbol{T}_2^{(t+1)} \boldsymbol{d}_2^{(t+1)} &= \bar{\boldsymbol{b}}^{(t+1)},
\f}
where the nodal displacements are separated into a vector of independent, unknown DOF's \f$\boldsymbol{d}_1\f$ and a vector of dependent DOF's \f$\boldsymbol{d}_2\f$. For a simulation using the displacement control approach without additional coupling equations, the matrix \f$\boldsymbol{T}_1^{(t+1)}\f$ vanishes and the matrix \f$\boldsymbol{T}_2^{(t+1)}\f$ is equal to the identity matrix with \f$\bar{\boldsymbol{b}}^{(t+1)}\f$ being the prescribed displacements. In the case of periodic boundary conditions (e.g. \f$d_k+d_j=0\f$) or in the case of coupling a 'hanging node' via displacement constraints (e.g. \f$d_l=0.5d_j+0.5d_k\f$), both transformation matrices do not vanish. These matrices depend on the load step indicated by \f$t\f$, since these conditions might change between successive load steps (e.g. \f$\bar{\boldsymbol{b}}\f$ is modified in each increment of a displacement controlled analysis). The assignment of a degree of freedom to either the dependent or the independent DOF's is sometimes arbitrary. For example in the case of the periodic boundary conditions, either \f$d_k\f$ or \f$d_j\f$ are independent, whereas the other is then the dependent DOF. Using a Gauss elimination procedure, which is equivalent to an inversion of the matrix \f$\boldsymbol T_2^{(t+1)}\f$, the vector of total displacements \f$\boldsymbol{d}\f$ can be expressed as a function of the independent DOF's \f$\boldsymbol{d}_1\f$ only:
\f{align*}{
 \boldsymbol{d}^{(t+1)}&=
\begin{bmatrix}
 \boldsymbol{d}_1^{(t+1)}\\
 \boldsymbol{d}_2^{(t+1)}
\end{bmatrix}
 =\begin{bmatrix}
                         \boldsymbol{I}\\
                         -\boldsymbol{C}_{con}
  \end{bmatrix}
 \boldsymbol{d}_1^{(t+1)}+
 \begin{bmatrix}
  \boldsymbol{0}\\
  \boldsymbol{b}^{(t+1)}
 \end{bmatrix},
\f}
where it is assumed that the vector \f$\boldsymbol d^{(t+1)}\f$ is reordered so that the independent parameters are at the top. Note the difference between the vector \f$\bar{\boldsymbol{b}}^{(t+1)}\f$  and \f$\boldsymbol{b}^{(t+1)}\f$. The matrix \f$\boldsymbol{C}_{con}\f$ relates the dependent DOF's with the independent DOF's. In the further derivation it is assumed that the matrix \f$C\f$ does not depend on the load step, since the increase of the loading is only obtained by modification of the right hand side \f$\boldsymbol{b}\f$.

Consequently, the variation of the nodal displacement vector is given by
\f{align*}{
\delta\boldsymbol{d}^{(t+1)}&
=\begin{bmatrix}
                        \boldsymbol{I}\\
                        -\boldsymbol{C}_{con}
 \end{bmatrix}
\delta\boldsymbol{d}_1^{(t+1)}
\f}

Using the element shape functions, the relation between nodal displacements \f$\boldsymbol{d}\f$ and strains \f$\boldsymbol{\varepsilon}\f$ can be expressed as
\f{align*}{
 \Delta\boldsymbol{\varepsilon}^{(t+1)}&=
 \boldsymbol{\varepsilon}^{(t+1)}- \boldsymbol{\varepsilon}^{(t)}\\
&=\boldsymbol{B}\left(
\begin{bmatrix}
\boldsymbol{I}\\
-\boldsymbol{C}_{con}
\end{bmatrix}
\left(\boldsymbol{d}_1^{(t+1)}-\boldsymbol{d}_1^{(t)}\right)
+
\begin{bmatrix}
\boldsymbol{0}\\
\boldsymbol{b}^{(t+1)}
\end{bmatrix}
-\begin{bmatrix}
\boldsymbol{0}\\
\boldsymbol{b}^{(t)}
\end{bmatrix}
\right)\\
\delta\boldsymbol{\varepsilon}^{(t+1)}
&=\boldsymbol{B}\left(
\begin{bmatrix}
\boldsymbol{I}\\
-\boldsymbol{C}_{con}
\end{bmatrix}
\delta\boldsymbol{d}_1^{(t+1)}
\right).
\f}
Setting
\f{align*}{
\Delta\boldsymbol{d}_1^{(t+1)}&=\left(\boldsymbol{d}_1^{(t+1)}-\boldsymbol{d}_1^{(t)}\right)
\f}
 and substituting of these equations into the linearized system of the virtual work equations gives
\f{align*}{
\left(
\begin{bmatrix}
\boldsymbol{I}\\
-\boldsymbol{C}_{con}
\end{bmatrix}
\Delta\boldsymbol{d}_1^{(t+1)}+
\begin{bmatrix}
\boldsymbol{0}\\
\boldsymbol{b}^{(t+1)}
\end{bmatrix}
-
\begin{bmatrix}
\boldsymbol{0}\\
\boldsymbol{b}^{(t)}
\end{bmatrix}
\right)^{T}
\boldsymbol{K}^T
\left(
\begin{bmatrix}
\boldsymbol{I}\\
-\boldsymbol{C}_{con}
\end{bmatrix}
\delta\boldsymbol{d}_1^{(t+1)}
\right)
= \left(\Delta \boldsymbol F_{ext}^{(t+1)}+\boldsymbol{F}_{ext}^{(t)}-\boldsymbol F_{int}^{(t)}\right)^T \begin{bmatrix}
                        \boldsymbol{I}\\
                        -\boldsymbol{C}_{con}
 \end{bmatrix}
\delta\boldsymbol{d}_1^{(t+1)}.
\f}
The integrated stiffness matrix and the internal and external force vectors are split into submatrices similar to \f$\boldsymbol d_1\f$ and \f$\boldsymbol d_2\f$:
\f{align*}{
\int_V \boldsymbol{B}^T\dfrac{\partial \sigma(\boldsymbol{d}^{(t)})}{\partial \varepsilon}\boldsymbol{B}\;dV&=\boldsymbol K=\begin{bmatrix}
  \boldsymbol K_{11} &  \boldsymbol K_{12}\\
   \boldsymbol K_{21} &   \boldsymbol K_{22}                                                                  \end{bmatrix}\\
\int_V\boldsymbol{B}^T\boldsymbol\sigma^{(t)}\;dV&=\boldsymbol F_{int}^{(t)}=
\begin{bmatrix}
\boldsymbol F_{1,int}^{(t)}\\
\boldsymbol F_{2,int}^{(t)}
\end{bmatrix}.
\f}
Canceling out \f$\delta\boldsymbol{d}_1^{(t+1)}\f$, transposing and splitting the external and internal forces into two parts yields:
\f{align*}{
\begin{bmatrix}
\boldsymbol{I}&
-\boldsymbol{C}_{con}^{T}
\end{bmatrix}
\begin{bmatrix}
  \boldsymbol K_{11} &  \boldsymbol K_{12}\\
   \boldsymbol K_{21} &   \boldsymbol K_{22}
\end{bmatrix}
\left(
\begin{bmatrix}
\boldsymbol{I}
-\boldsymbol{C}_{con}
\end{bmatrix}
\Delta\boldsymbol{d}_1^{(t+1)}+
\begin{bmatrix}
\boldsymbol{0}\\
\boldsymbol{b}^{(t+1)}
\end{bmatrix}
-\begin{bmatrix}
\boldsymbol{0}\\
\boldsymbol{b}^{(t)}
\end{bmatrix}
\right)
=\begin{bmatrix}
                        \boldsymbol{I}&
                        -\boldsymbol{C}^T_{con}
 \end{bmatrix}
\begin{bmatrix}
\Delta \boldsymbol F_{1,ext}^{(t+1)}+\boldsymbol F_{1,ext}^{(t)}-\boldsymbol F_{1,int}^{(t)}\\
\Delta \boldsymbol F_{2,ext}^{(t+1)}+\boldsymbol F_{2,ext}^{(t)}-\boldsymbol F_{2,int}^{(t)} 
\end{bmatrix}.
\f}
Expanding this equation finally gives:
\f{align*}{
\boldsymbol{K}^{mod}\Delta\boldsymbol{d}_1^{(t+1)}=\left(\boldsymbol{K}_{12}-\boldsymbol{C}^T_{con}\boldsymbol{K}_{22}\right)\left(\boldsymbol{b}^{(t)}-\boldsymbol{b}^{(t+1)}\right)+\Delta \boldsymbol F_{1,ext}^{(t+1)}+\boldsymbol F_{1,ext}^{(t)}\\
-\boldsymbol F_{1,int}^{(t)}
-\boldsymbol{C}^T_{con}\left(\Delta \boldsymbol F_2^{(t+1)}+\boldsymbol F_{2,ext}^{(t)}-\boldsymbol F_{2,int}^{(t)}
\right)
\f}
with the modified stiffness matrix
\f{align*}{
\boldsymbol K^{mod}&= \boldsymbol{K}_{11}-\boldsymbol{C}^T_{con}\boldsymbol{K}_{21}-\boldsymbol{K}_{12}\boldsymbol{C}_{con}+\boldsymbol{C}_{con}^T\boldsymbol{K}_{22}\boldsymbol{C}_{con}.
\f}

# Implementation

## Basics 

- The class `NuTo::Constraint::Constraints` has a container of `NuTo::Constaint::Equation`. 
- Each `NuTo::Constraint::Equation` has a time dependent right-hand-side (`std::function<double(double time)>`) and a container of `NuTo::Constraint::Term`. 
- Each term saves a reference to a `NuTo::NodeBase`, a component and a coefficient.

That way, `NuTo::Constraint::Constraints` can build the matrix \f$\boldsymbol T\f$ and the rhs \f$\boldsymbol b\f$ of

\f[
\boldsymbol T \boldsymbol u = \boldsymbol b(t)
\f]

The method `NuTo::Assembler::BuildGlobalDof()` performs a sparse (pivot?) Gauss elimination and identifies the dependent dofs \f$\boldsymbol u_2\f$, performs renumbering, and so on... After that, the following members are valid:

- `BlockSparseMatrix mConstraintMatrix` - corresponds to \f$\boldsymbol T_2^{-1} \boldsymbol T_1\f$
- `BlockSparseMatrix mConstraintRhsMapping` - corresponds to \f$\boldsymbol T_2^{-1}\f$
- `BlockFullVector mConstraintRhs` - corresponds to \f$\boldsymbol T_2^{-1} \boldsymbol b(t)\f$

`mConstraintRhs` is updated to a new global time via `NuTo::ConstraintUpdateRhs(time)` which internally calculates \f$\boldsymbol b(t)\f$ and applies \f$\boldsymbol T_2^{-1}\f$ (which is equal, but way faster, than performing the Gauss elimination again (\f$O(n^3)\f$?) - at every time step.)


## Usage

1) Create as many instances of `NuTo::Constraint::Equation` as you need
    - manually by instantiating an object of that class (not recommended)
    - using `mechanics/constraints/ConstraintCompanion.h`
        - Just have a look at this class, very intuitive, always returns an `NuTo::Constraint::Equation` or a `std::vector` of those

2) Store the equation
    - call `NuTo::Constraint::Constraints::Add(...)`
    - you can access this class via `NuTo::StructureBase::Constraints()`

