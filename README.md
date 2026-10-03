# IRP finite-volume-element scheme for polydisperse flow and segregation

This Github repository contains the source files of a finite volume element solver written in MATLAB designed to approximate the coupled  transport-flow problem:

$$
	\begin{aligned}
			\partial_t  \phi_{ l} +\nabla \cdot \big(\phi_{ l} \boldsymbol{u}+  f_{ l}(\Phi)\boldsymbol{k}\big)&=0,\quad  
		f_l ( \Phi) = \phi_l v_l ( \Phi), \quad l=1,\dots,N; \\
			   -\nabla \cdot \bigl(\mu(\phi)\boldsymbol{\varepsilon}(\boldsymbol{u}) \bigr)+ \nabla p &= \boldsymbol{g}(\phi),
		\quad \boldsymbol{\varepsilon} ( \boldsymbol{u} ) \coloneqq \frac12 ( \nabla \boldsymbol{u} 
		+ (\nabla \boldsymbol{u})^{\mathrm{T}}); 
		\\
			 \quad \nabla \cdot \boldsymbol{u} &= 0, 
	\end{aligned}
 $$ 

where  $`\Phi \coloneqq ( \phi_1, \dots, \phi_N)^{\mathrm{T}}`$, $`\phi \coloneqq \phi_1 + \cdots + \phi_N`$,  and is posed on a bounded domain $`\Omega \subset \mathbb{R}^2`$ for $`t >0`$. The system can be  understood as a transport-flow model for a  two-phase mixture consisting of a disperse phase with $`N`$ species of particles or droplets moving  in a viscous continuous phase. In particular it models  a polydisperse suspension    of solid particles of $`N`$ species with diameters $`d_1 \geq \dots \geq d_N`$ and densities $`\rho_1, \dots, \rho_N`$  dispersed   in a viscous fluid. The  unknowns are the volume fractions $`\phi_l=\phi_l (\boldsymbol{x}, t)`$ of each particle species $`l`$ (having diameter $`d_l`$ and density $`\rho_l`$), which depend on  spatial position $`\boldsymbol{x}`$ and time $`t`$, the volume-averaged mixture velocity $`\boldsymbol{u}=\boldsymbol{u}(\boldsymbol{x},t)`$,  and the pressure $`p=p(\boldsymbol{x},t)`$. Here  $`\boldsymbol{k}`$ denotes the downward-pointing unit vector   while $v_1(\Phi), \dots, v_N(\Phi)$ are  prescribed velocity functions 
 that  describe the relative motion of each particle species with respect to the   mixture. The physically relevant set of states $`\Phi`$ is 
 
$$
\begin{aligned}
D := \{  (  \phi_1, \dots, \phi_N)^{\mathrm{T}} \in \mathbb{R}^N  :  
	\phi_1 \geq 0, \dots,\phi_N \geq 0,  \phi := \phi_1 + \cdots + \phi_N \leq \phi_{\max}  \},
	\end{aligned}
	$$

where $\phi_{\max}$ denotes a given maximum total particle volume fraction. 

This FVE  scheme is second-order accurate in both space and time and produce numerical solutions satisfying the invariant region preservation (IRP) property, that is, if $`\Phi_0 (\boldsymbol{x}) \in \mathcal{D}`$ for all $`\boldsymbol{x} \in \Omega`$, then the numerical solution assumes 
  values in $`\mathcal{D}`$ for all times. The Stokes problem is discretised with nonconforming Crouzeix–Raviart/P0
elements on a triangulation, and the $`N`$
concentration equations with a finite volume scheme on the dual diamond
mesh: LLF flux, MUSCL reconstruction with scaling limiters, and SSPRK2 time
stepping (Algorithm 3.1).

We kindly ask you to to acknowledge the use of this software by citing the (current) paper:

> J. Barajas-Calonge, R. Bürger, P. Mulet, L. M. Villada,
> *An invariant-region-preserving finite-volume-element scheme for a
> polydisperse model of flow and segregation in general domains* (2026).

Related papers presenting high-order IRP schemes for the one-dimensional version of the model and its two-dimensional version on Cartesian grids:

- J. Barajas-Calonge, R. Bürger, P. Mulet and L.M. Villada. **Invariant-region-preserving WENO schemes for one-dimensional multispecies kinematic flow models**, 
*J. Comput. Phys.* 537 (2025), article 114081. [**link**](https://www.sciencedirect.com/science/article/abs/pii/S002199912500364X)
- J. Barajas-Calonge, R. Bürger, P. Mulet and L.M. Villada. **A second-order invariant-region-preserving scheme for a transport-flow model of polydisperse sedimentation**,
*J. Sci. Comput.* 108 (2026), article 29, [**link**](https://link.springer.com/article/10.1007/s10915-026-03356-y)


## Contents

```
ex1_accuracy.m           Example 1  manufactured solution 
ex2_rotated_vessels.m    Example 2  FVE vs IRP-FVE in rotated vessels 
ex3_roof.m               Example 3  roof-shaped vessel 
ex4_reflux_classifier.m  Example 4  reflux classifier
ex5_inflow_outflow.m     Example 5  tank with inflow and outflow 
plot_result.m            generic plot of a stored snapshot
src/                     solver
Meshes/                  meshes (.mat)
Results/                 snapshots written by the examples (.mat)
```

`src/`:

| file | purpose |
|---|---|
| `default_model.m` | model parameters of Section 4 (MLB model, viscosity, buoyancy) |
| `build_mesh.m`, `load_mesh.m` | primal triangulation → dual mesh; reading the mesh files |
| `stokes_setup.m`, `stokes_solve.m` | CR/P0 Stokes solver, problem (2.2) |
| `fv_setup.m`, `fv_residual.m` | dual-mesh FV scheme: MUSCL gradient, slope and scaling limiters, LLF flux |
| `solve_fve.m` | coupled time loop (Algorithm 3.1), CFL condition (3.12), snapshots |
| `cell_average.m` | cell averages of initial data |
| `manufactured_ex1.m` | exact solution and forcing terms of Example 1 |
| `inflow_outflow_bc.m` | boundary data of Example 5 |

## Usage

Run the example scripts from the repository root, e.g.

```matlab
>> ex3_roof
```

The scripts only compute; tables and figures are obtained from the stored
snapshots. Snapshots are named `Results/<example>_..._t<time>.mat` and contain

| variable | content |
|---|---|
| `u` | cell averages on the dual mesh, NC x N |
| `qx`, `qy` | Crouzeix–Raviart velocity (mean value on each primal edge) |
| `p` | pressure, one value per triangle |
| `t` | time |
| `minphi`, `maxphi` | `min_{K,n} phi_{l,K}^n` (1 x N) and `max_{K,n} phi_K^n` over all steps up to `t` |
| `g` | dual mesh: vertices `V`, cells `T` (4 x NC, `T(4,:) = -1` on boundary triangles), areas `A`, centres `xC`, `yC`, faces `F`, `N`, `M`, `L` |
| `gt` | primal mesh: vertices `V`, triangles `T`, barycentres `B`, areas `A`, edges `F`, `N`, `M`, `L`, and `elem2dof` (edges of each triangle) |

One dual control volume is associated with each edge of the triangulation,
with the same index, so `u(e,:)`, `qx(e)` and `qy(e)` all refer to edge `e`.

A snapshot can be plotted with

```matlab
>> plot_result('Results/ex3_k4_irp_t3.mat', 'phi', 'view', [90 90])
>> plot_result('Results/ex4a_t4.mat', 'umag', 'rotate', 90, 'view', [90 90], 'quiver', 10)
>> plot_result('Results/ex5_t4p25.mat', 'phi1', 'export', 'ex5_phi1.eps')
```

(fields `phi1`, ..., `phiN`, `phi`, `ux`, `uy`, `umag`, `p`; see
`help plot_result`).

Mesh files: `square<n>.mat` (Example 1, n = 4, ..., 128),
`inclined_<theta>_<k>.mat` (Example 2), `roof_<k>.mat` (Example 3),
`refluxclas_<k>.mat` (Example 4) and `tank_160x30.mat` (Example 5). A mesh file needs at least the triangulation
`Vtri` (2 x nv) and `Ttri` (3 x nt); see `load_mesh.m`.

