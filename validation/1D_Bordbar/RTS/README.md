# Radiative Transfer Simulator

![](https://github.com/gusirosx/RTS/blob/master/RTS.jpg)

RTS or Radiative Transfer Simulator is the specialized radiation code initially developed and designed by [Gustavo Silva Rodrigues](https://github.com/gusirosx). This code was written to solve the stationary form of the radiative transport equation in Cartesian meshes in uni, bi and three-dimensional computational domains. In addition, RTS is being developed considering the necessary flexibility that every radiation code must have. The purpose of developing this code is to provide a modular tool that can be coupled with the vast majority of CFD algorithms. It should be noted that anyone can use, copy, distribute, study, change and improve it according to General Public License ([GNU GPL v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)). 

Among the main characteristics of the RTS, it is worth mentioning:

* Implemented in Fortran 90;
* It has a modular structure;
* Supports 1D, 2D or 3D orthogonal Cartesian meshes;
* You can solve problems using non-uniform Cartesian meshes;
* You can solve the RTE with the discretizing schemesupwind and CDS;
* They have directional solvers P1, DOM and FAM
* Uses the quadrature sets Sn, Tn and Qn in the Discrete Ordinates Method;
* Can solve problems with an emitting, absorbent and dispersive medium;
* Support for WSGG (Weighted Sum of Gray Gases Model) absorption models in gray and non-gray formulation;
* Can solve problems with isotropic and anisotropic scattering;
* Support Mie and Rayleigh scattering functions;
* Able to solve problems with non-homogeneous medium;
* Able to solve problems with black walls or diffuse reflections;
* Able to solve problems in radiative equilibrium;
* Does not depend on any external library, being a pure Fortran program.
