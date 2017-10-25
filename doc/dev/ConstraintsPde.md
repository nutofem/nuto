@page ConstraintsPde Handling of constraints in NuToPde

\f$\newcommand\vec[1]{\boldsymbol{#1}}\f$
# Theory

From the principle of virtual displacements, it follows for a body in equilibrium:
\f[
\delta W_{int}^{(t+1)} = \delta W_{ext}^{(t+1)},
\f]
where the variation of the internal energy \f$\delta W_{int}^{(t+1)}\f$ and the
variation of the external energy \f$\delta W_{ext}^{(t+1)}\f$ is given by
\f{align*}{
    \delta W_{int}^{(t+1)}
    &=
    \int_V\vec{\sigma}^T(\vec{d}^{(t+1)})\delta\vec{\varepsilon}^{(t+1)}\;dV\\
    &=
    \int_V\left(\vec{\sigma}(\vec{d}^{(t)})+\Delta \vec\sigma\left(\vec{d}^{(t)}+\Delta \vec{d}^{(t+1)}\right)\right)^T\delta\vec\varepsilon^{(t+1)}\;dV
    \\
    \delta W_{ext}^{(t+1)}
    &=
    \left(\vec F_{ext}^{(t)}+\Delta \vec{F}_{ext}^{(t+1)}\right)^T \delta \vec{d}^{(t+1)}.
\f}

Linearization of the system of equations gives
\f[
    \int_V
    \left(
        \dfrac{\partial \sigma(\vec{d}^{(t)})}{\partial \varepsilon}\Delta\vec\varepsilon^{(t+1)}
    \right)^T
    \delta\vec\varepsilon^{(t+1)}\;dV
    =
    \left(
        \Delta \vec F_{ext}^{(t+1)}
    \right)^T
    \delta \vec{d}^{(t+1)}
    -
    \int_V
    \left( \vec\sigma^{(t)} \right)^T
    \delta\vec\varepsilon^{(t+1)}\;dV
    +
    \left( \vec{F}_{ext}^{(t)} \right)^T
    \delta\vec d^{(t+1)},
\f]
where the last two summands on the right hand side vanish, if the previous step
is in equilibrium.  In the general formulation, the vector \f$\vec{d}\f$
includes all the nodal displacements. However, if kinematic boundary conditions
have to be applied, the variations of the displacement state have to be
kinematically compatible. In general, the linear constraint equations can be
written as
\f[
\vec{T}_1^{(t+1)} \vec{d}_1^{(t+1)}+\vec{T}_2^{(t+1)} \vec{d}_2^{(t+1)} = \vec{b}^{(t+1)},
\f]
where the nodal displacements are separated into a vector of independent,
unknown DOFs \f$\vec{d}_1\f$ and a vector of dependent DOFs
\f$\vec{d}_2\f$. 

We restrict the way, constraint equations can be build.
 - Each constraint equation must only contain exactly one dependent dof. If you need to 
   define a constraint equation with more than one dependent dof, you have to perform an 
   elimination manually.
 - The coefficient in front of the dependent dof term has to be `1`.

With the correct DOF numbering, this ensures that \f$\vec{T}_2\f$ is an identity matrix 
and we can trivially eliminate the dependent DOFs \f$\vec{d}_2\f$.

\f[
    \vec{d}^{(t+1)}
    =
    \begin{bmatrix}
        \vec{d}_1^{(t+1)}\\
        \vec{d}_2^{(t+1)}
    \end{bmatrix}
    =
    \begin{bmatrix}
        \vec{I}\\
        -\vec{T}_{1}
    \end{bmatrix}
    \vec{d}_1^{(t+1)}
    +
    \begin{bmatrix}
        \vec{0}\\
        \vec{b}^{(t+1)}
    \end{bmatrix},
\f]

The matrix \f$\vec{T}_1\f$ is the so called _constraint matrix_ and will now be called 
\f$\vec{C}_{con} \f$.

\f[
    \vec{d}^{(t+1)}
    =
    \begin{bmatrix}
        \vec{d}_1^{(t+1)}\\
        \vec{d}_2^{(t+1)}
    \end{bmatrix}
    =
    \begin{bmatrix}
        \vec{I}\\
        -\vec{C}_{con}
    \end{bmatrix}
    \vec{d}_1^{(t+1)}
    +
    \begin{bmatrix}
        \vec{0}\\
        \vec{b}^{(t+1)}
    \end{bmatrix}
\f]

The matrix \f$\vec{C}_{con} \f$ relates the dependent DOFs with the
independent DOFs. In the further derivation it is assumed that the matrix
\f$\vec{C}_{con}\f$ does not depend on the load step, since the increase of the loading is
only obtained by modification of the right hand side \f$\vec{b}\f$.

Consequently, the variation of the nodal displacement vector is given by
\f[
    \delta\vec{d}^{(t+1)}
    =
    \begin{bmatrix}
        \vec{I}\\
        -\vec{C}_{con}
    \end{bmatrix}
    \delta\vec{d}_1^{(t+1)}.
\f]

Using the element shape functions, the relation between nodal displacements
\f$\vec{d}\f$ and strains \f$\vec{\varepsilon}\f$ can be
expressed as
\f{align*}{
    \Delta\vec{\varepsilon}^{(t+1)}
    &=
    \vec{\varepsilon}^{(t+1)}- \vec{\varepsilon}^{(t)}\\
    &=
    \vec{B}\left(
        \begin{bmatrix}
            \vec{I}\\
            -\vec{C}_{con}
        \end{bmatrix}
        \left(
            \vec{d}_1^{(t+1)}-\vec{d}_1^{(t)}
        \right)
        +
        \begin{bmatrix}
            \vec{0}\\
            \vec{b}^{(t+1)}
        \end{bmatrix}
        -
        \begin{bmatrix}
            \vec{0}\\
            \vec{b}^{(t)}
        \end{bmatrix}
        \right)\\
    \delta\vec{\varepsilon}^{(t+1)}
    &=
    \vec{B}\left(
        \begin{bmatrix}
            \vec{I}\\
            -\vec{C}_{con}
        \end{bmatrix}
        \delta\vec{d}_1^{(t+1)}
    \right).
\f}

Setting
\f[
\Delta\vec{d}_1^{(t+1)}=\left(\vec{d}_1^{(t+1)}-\vec{d}_1^{(t)}\right)
\f]
and substituting of these equations into the linearized system of the virtual work equations gives
\f[
    \left(
        \begin{bmatrix}
            \vec{I}\\
            -\vec{C}_{con}
        \end{bmatrix}
        \Delta\vec{d}_1^{(t+1)}+
        \begin{bmatrix}
            \vec{0}\\
            \vec{b}^{(t+1)}
        \end{bmatrix}
        -
        \begin{bmatrix}
            \vec{0}\\
            \vec{b}^{(t)}
        \end{bmatrix}
    \right)^{T}
    \vec{K}^T
    \left(
        \begin{bmatrix}
            \vec{I}\\
            -\vec{C}_{con}
        \end{bmatrix}
        \delta\vec{d}_1^{(t+1)}
    \right)
    =
    \left(
        \Delta \vec F_{ext}^{(t+1)}
        +
        \vec{F}_{ext}^{(t)}
        -
        \vec F_{int}^{(t)}
    \right)^T
    \begin{bmatrix}
        \vec{I}\\
        -\vec{C}_{con}
    \end{bmatrix}
    \delta\vec{d}_1^{(t+1)}.
\f]
The integrated stiffness matrix and the internal and external force vectors are
split into submatrices similar to \f$\vec d_1\f$ and \f$\vec
d_2\f$:
\f{align*}{
    \int_V \vec{B}^T\dfrac{\partial \sigma(\vec{d}^{(t)})}{\partial \varepsilon}\vec{B}\;dV
    &=
    \vec K
    =
    \begin{bmatrix}
        \vec K_{11} &  \vec K_{12}\\
        \vec K_{21} &  \vec K_{22}
    \end{bmatrix}\\
    \int_V\vec{B}^T\vec\sigma^{(t)}\;dV&=\vec F_{int}^{(t)}
    =
    \begin{bmatrix}
        \vec F_{1,int}^{(t)}\\
        \vec F_{2,int}^{(t)}
    \end{bmatrix}.
\f}
Canceling out \f$\delta\vec{d}_1^{(t+1)}\f$, transposing and splitting the external and internal forces into two parts yields:
\f[
    \begin{bmatrix}
        \vec{I}&
        -\vec{C}_{con}^{T}
    \end{bmatrix}
    \begin{bmatrix}
        \vec K_{11} &  \vec K_{12}\\
        \vec K_{21} &  \vec K_{22}
    \end{bmatrix}
    \left(
        \begin{bmatrix}
            \vec{I}
            -\vec{C}_{con}
        \end{bmatrix}
        \Delta\vec{d}_1^{(t+1)}
        +
        \begin{bmatrix}
            \vec{0}\\
            \vec{b}^{(t+1)}
        \end{bmatrix}
        -
        \begin{bmatrix}
            \vec{0}\\
            \vec{b}^{(t)}
        \end{bmatrix}
    \right)
    =
    \begin{bmatrix}
        \vec{I}&
        -\vec{C}^T_{con}
    \end{bmatrix}
    \begin{bmatrix}
        \Delta \vec F_{1,ext}^{(t+1)}+\vec F_{1,ext}^{(t)}-\vec F_{1,int}^{(t)}\\
        \Delta \vec F_{2,ext}^{(t+1)}+\vec F_{2,ext}^{(t)}-\vec F_{2,int}^{(t)} 
    \end{bmatrix}.
\f]
Expanding this equation finally gives:
\f{align*}{
    \vec{K}^{mod}\Delta\vec{d}_1^{(t+1)}
    =
    \left(
        \vec{K}_{12}-\vec{C}^T_{con}\vec{K}_{22}
    \right)
    \left(
        \vec{b}^{(t)}-\vec{b}^{(t+1)}\right)+\Delta \vec F_{1,ext}^{(t+1)}+\vec F_{1,ext}^{(t)}\\
        -\vec F_{1,int}^{(t)} -\vec{C}^T_{con}\left(\Delta \vec F_2^{(t+1)}+\vec F_{2,ext}^{(t)}-\vec F_{2,int}^{(t)}
    \right)
\f}
with the modified stiffness matrix
\f[
\vec K^{mod} = \vec{K}_{11}-\vec{C}^T_{con}\vec{K}_{21}-\vec{K}_{12}\vec{C}_{con}+\vec{C}_{con}^T\vec{K}_{22}\vec{C}_{con}.
\f]

# Implementation

## Basics 

- The class `NuTo::ConstraintPde::Constraints` has a container of `NuTo::Constaint::Equation`. 
- Each `NuTo::ConstraintPde::Equation` has a time dependent right-hand-side (`std::function<double(double time)>`) and a container of `NuTo::ConstraintPde::Term`. 
- Each term saves a reference to a `NuTo::NodeSimple`, a component and a coefficient.
- The first term of each equation is the dependent term. Its coefficient is always `1`.
- The method `NuTo::ConstraintPde::Constraints::BuildConstraintMatrix(...)` returns a sparse matrix containing \f$\vec{C}_{con}\f$
- The method `NuTo::ConstraintPde::Constraints::GetRhs()` returns a the right hand side of the equations \f$\vec{b}\f$ 
- The constraints do not include a numbering of the dofs in `NuTo::NodeSimple`. This is done via `NuTo::DofNumbering::Build()` and has to be performed before you try to build \f$\vec{C}_{con}\f$

## Usage

1) Create as many instances of `NuTo::ConstraintPde::Equation` as you need
    - manually by instantiating an object of that class (not recommended)
    - using `mechanics/constraintsPde/ConstraintCompanion.h`
        - Just have a look at this class, very intuitive, always returns an `NuTo::ConstraintPde::Equation` or a `std::vector` of those

2) Store the equation
    - call `NuTo::ConstraintPde::Constraints::Add(...)`

3) Perform some kind of dof numbering
    - manual
    - using `NuTo::DofNumbering::Build()`

