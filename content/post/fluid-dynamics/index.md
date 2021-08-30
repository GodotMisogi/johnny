---
date: "2021-05-14T00:00:00+08:00"
title: A Flurry of Fluid Dynamics
draft: true
description: A quick introduction to fluid dynamics.
image: 
slug: fluid-dynamics
categories: 
    - Aerodynamics
    - Computational Fluid Dynamics
    - Programming
    - Mathematics
---

#### References

1. Batchelor, G. K. _An Introduction to Fluid Dynamics._ Cambridge University Press, 1992.
2. Drela, Mark. _Flight Vehicle Aerodynamics_. The MIT Press, 2014.
3. Katz, Joseph, Plotkin, Allen. _Low-Speed Aerodynamics - 2nd Edition._ Cambridge University Press, 2001.
4. Feynman, Richard P., Leighton, Robert B., Sands, Matthew. _The Feynman Lectures on Physics, Vol. 2: Mainly Electromagnetism and Matter_. Basic Books, 2011.
5. Balakrishnan, V. _Selected Topics in Mathematical Physics_. NPTEL, 2014.
 
## Definitions

The **circulation** $\Gamma$ of a vector field $\vec f(x^i)$ about a closed contour $C$ parametrised by a vector $d x_i$:

$$ \Gamma \equiv \oint_C \vec f(x^i)\ dx_i $$

Note that the dimension of the vector space is irrelevant.

This is related to **Stokes' theorem**, which indicates the 'total' curl of the vector field for an open surface $\partial \mathcal S$:

$$ \oint_C \vec f \cdot d\vec l = \iint_{\partial \mathcal S} \left(\nabla \times \vec f\right)\cdot \vec n\ d\mathcal S $$\

Similarly, the total 'emittance' of $\vec f$ over a closed surface $S$ is related to the 'total' divergence of a vector field from a volume $\mathcal V$ via **Gauss' theorem**:

$$ \iint_{\mathcal S} \vec f \cdot d\vec {\mathcal S} = \iiint_{\mathcal V}\nabla \cdot \vec f \ d\mathcal V $$

_Note:_ Gauss' and Stokes' theorems are properties of integrals, and are not particularly related to vector spaces as such; they are independent of tensor dimensionality.

The **material/Lagrangian/etc... derivative** $D/Dt$ is an operator expressing the chain rule on any tensor field that signifies the time rate of change of the tensor if it is dependent on space-time macroscopic velocity field variations. In a flat spacetime:

$$ \frac{D}{Dt} \equiv \frac{\partial}{\partial t} + v^i\frac{\partial}{\partial x_i} $$

This can be defined more generally on a smooth manifold by upgrading the spatial partial derivatives to covariant derivatives.

The conservation of an _extensive property_ $N$ as it travels through a defined **control volume** $\mathcal V$, in which the conservation is described in terms of an _intensive property_ $\eta$ (such as the density or any equivalent macroscopic property corresponding to a continuum formulation), is dictated by the **Reynolds' transport theorem** in the following temporal formulation:

$$ \frac{DN}{Dt} = \iiint_{\mathcal V} \frac{\partial \eta}{\partial t}\ d\mathcal{V} + \iint_{\mathcal S} \eta \vec V \cdot d\vec{\mathcal S} $$

The surface integral can usually be converted into a differential form over the volume via Gauss' theorem.

#### Fluid  Dynamics

A **fluid** is defined as a continuum of particles such that its microscopic properties average into macroscopic properties that approximate its behaviour as a single element, a behaviour characterised by negligibly resistive non-retention of shape when deformed.

The definitions of circulation and the material derivative are sufficient to derive **Kelvin's circulation theorem** for velocity fields $\vec V$ along a fluid curve $C$, defined as a curve which always passes through the same fluid elements:

$$ \frac{D\Gamma}{Dt} = \frac{D}{Dt}\oint_C \vec V \cdot d\vec l = \oint_C \frac{D\vec V}{Dt} \cdot d\vec l + \oint_C \vec V \cdot \frac{D}{dt} d\vec l$$

$$ \implies \frac{D\Gamma}{Dt} = \oint_C \vec a \cdot d\vec l, \quad \because \frac{D\vec V}{Dt} = \vec a,\quad \oint_C \vec V \cdot d\vec V = 0 $$

The differential form of the **continuity equation** representing the conservation of mass $M$, with a macroscopic property called the **density** $\rho$, via the Reynolds' transport theorem is:

$$ \frac{DM}{Dt} = \frac{\partial \rho}{\partial t} + \nabla\cdot (\rho \vec V) $$

For a control volume with no time-varying source or sink of mass, $dM/dt = 0$, which gives the usual continuity equation analysed in fluid mechanics. An **incompressible flow** is defined as a flow in which the changes in density with respect to space and time are negligible.

The differential form of the **incompressible Navier-Stokes equations** in an external time-independent gravitational field $\vec g$ with non-conservative body forces $\vec f$ is:
    
$$ \nabla \cdot \vec V = 0 $$

$$ \rho \frac{D\vec V}{Dt} = -\nabla p + \rho \vec g + \nabla \cdot \hat \tau + \rho\vec f$$

The energy equation reduces to Bernoulli's equation (mentioned later), which can be derived from the momentum equation above, hence it need not be considered. 

For a Newtonian fluid, the stresses are considered to be linearly proportional to the velocity gradient near a wall, with a constant $\mu$ called **viscosity**:

$$ \tau = \mu \frac{\partial u}{\partial x} $$ 

Hence the stress tensor is (explicitly showing the Kronecker product $\otimes$):

$$\hat \tau = \mu\left[\nabla \otimes \vec V +  \left(\nabla \otimes \vec V\right)^T  - \frac{2\nu}{3} \hat I\right]$$

where $\nu \equiv \nabla \cdot \vec V$, called the **dilatation strain rate**.

Non-dimensionalisation schemes are important in fluid dynamics because they remove the dependence of the equations on intrinsic fluid properties, such as $\rho,\ \mu$ and reference length scales such as $L$ and $|\vec V|$ for flow over objects.

If the flow is inviscid, the time derivatives and the stress terms drop out, giving the **incompressible Euler equations**:

$$ \nabla \cdot \vec V = 0 $$

$$ \rho(\vec V \cdot \nabla \vec V) = -\nabla p + \rho\vec g + \vec f$$

Assume that the external gravitational field is conservative, and hence is expressible as a scalar potential $\vec g = -\nabla \xi$, as is usually done in physics. Assuming no non-conservative body forces, Kelvin's circulation theorem provides the following statement of angular momentum conservation:

$$ \frac{D\Gamma}{Dt} = \oint_C \nabla\left(\frac{-p}{\rho}\right) \cdot d\vec l + \frac{1}{\rho}\oint_C \vec g \cdot d\vec l = 0 $$

The **vorticity** $\omega$ is defined as the curl of the velocity field $\vec V$, specifically: $ \omega \equiv \nabla \times \vec V $. Assume the flow is irrotational, i.e. $\omega = 0$, and there are no body forces, because the assumptions aren't enough. Using Stokes' theorem, $\Gamma = 0$. This describes a _conservative_ vector field, which can be expressed in the form of a scalar potential $\vec V = \nabla \phi$, when inserted into the continuity equation gives **Laplace's equation**:

$$ \nabla^2 \phi = 0 $$

As the real and imaginary parts of an analytic function of a complex variable $z$ independently satisfy Laplace's equation, another Laplace equation can be constructed by considering $\phi$ as the real part and using the Cauchy-Riemann conditions in Cartesian coordinates:

$$ \vec V = (u(x,y), v(x,y)), \quad u = \frac{\partial \phi}{\partial x} = \frac{\partial \psi}{\partial y},\quad v = \frac{\partial \phi}{\partial y} = -\frac{\partial \psi}{\partial x} $$

where $\psi$ is called a **stream function**, mainly used in two-dimensional flows. This leads to the construction of a **complex potential** $\Psi(z) = \phi + i \psi $. The derivative of this with respect to $z$ gives you the **complex velocity** $u - iv$. The stream function can also be constructed by the following definition of a **streamline** in three dimensions, a line which satsfies the equations:

$$ \vec V \times d\vec l = 0 $$ 

In two dimensions, say the $x$-$y$ plane, which physically relates the changes in velocity with respect to changes in space:

$$ \frac{dy}{dx} = \frac{v(x,y)}{u(x,y)} \iff u\ dy - v\ dx = 0$$

The second expression is an exact differential $d\psi =0$ if and only if the previous conditions for $\psi$ are satisfied. Note that the contours generated by $\phi$ and $\psi$ are orthogonal.

**Vortex lines** can be defined analogously to the definition of the streamline:

$$ \omega \times d\vec l = 0 $$

**Non-dimensionalisation** is a key concept in which the equations of motion are multiplied/divided by reference scales corresponding to certain physical quantities in external flows. In general, these are the length, speed, time, pressure and body forces.

$$ x^{i*} = \frac{x^i}{L}, \ u^{i*} = \frac{u^i}{V}, \ t^* = \frac{t}{T}, \ p^* = \frac{p}{p_0}, \ f^{i*} = \frac{f^i}{f_0} $$

The corresponding derivatives in non-dimensional terms can be obtained via the chain rule, which results in the equations of motion as follows:

$$ \frac{V}{L}\left(\frac{\partial u^{i*}}{\partial x^{i*}}\right) = 0 $$

#### Aside: Mathematical Generality

Electromagnetics and fluid dynamics are usually about almost exclusively solving partial differential equations, either analytically or numerically. The Laplace equation also governs the basic laws of electrostatics, which can be derived from Maxwell's equations. So let's look at some basic solutions first for the **Poisson equation** with a source $\rho$, a more general form of Laplace's equation first in $N$ dimensions:

$$ \nabla_N^2 \phi(\vec x) = \rho(\vec x) $$

This can be solved via the method of Green functions:

$$
\begin{aligned}
\nabla_N^2 G(\vec x, \vec x_0) & = \delta(\vec x - \vec x_0) \\\\ 
\int_{-\infty}^{\infty} \left[\nabla_N^2 G(\vec x, \vec x_0)\right]  e^{i\vec k \cdot \vec x}\ d\vec x & = \int_{-\infty}^{\infty} \delta(\vec x - \vec x_0) e^{i\vec k \cdot \vec x}\ d\vec x \\\\ 
|\vec k|^2 G(\vec k, \vec x_0) & = e^{i \vec k \cdot \vec x_0} \\\\ 
\implies G(\vec x, \vec x_0) & = \frac{1}{(2\pi)^N}\int_{-\infty}^{\infty} \frac{e^{-i \vec k \cdot (\vec x -\vec x_0)}}{|\vec k|^2} \ d\vec k 
\end{aligned}
$$

Laplace's equation is spherically symmetric in $N$ dimensions (evident from the $|\vec k|^2$ term), so expressing the integral in a $N$-spherical polar coordinate system will be more convenient. This is singular with multiple poles at $\vec k = \vec 0$, so the contour must be deformed in the complex plane.

$$ G(\vec x, \vec x_0) = \frac{1}{(2\pi)^N} \int_0^{2\pi} \int_{0}^{\infty} \frac{e^{-ik_r r\cos\theta}}{|r|^2}r\ dr d\theta $$

#### Potential Flows

So a majority of basic (and unrealistic!) fluid flows can be generated by solving Laplace's equation for $\psi$ (the streamlines) using certain boundary conditions and declarations of $u, v$, which is somewhat of an inverse problem. The incompressible Euler equation is then used to evaluate the pressure. Sometimes coordinate transformations would be required to integrate the PDEs more easily, but they might introduce scale factors. Note that the cases analysed are 2D, but can be generalised readily. As von Neumann said, analysis under so many assumptions makes this turns this exercise into the study of 'dry water'.

Laplace's equation is linear, so any superposition of different solutions will give you a new function describing a more general flow. A common point of analysis in potential flows is the **stagnation point**, at which $u(x,y) = v(x,y) = 0$. This is equivalent to finding the extrema of the solutions to the potential or stream functions, and the minima correspond to stagnation points. The streamline corresponding to the stagnation point is where the velocity vanishes, called a **dividing streamline**. It can be considered as a boundary with no normal velocity component of a specific shape, which corresponds to the solution of a flow around this boundary. This is a Neumann boundary condition:

$$ \nabla \phi \cdot \hat n = 0 $$

Usually in aerodynamic flows, the **pressure coefficient** is an important quantity to measure local pressure $p$ at some point compared to the freestream pressure $p_{\infty}$, using the Bernoulli equation in an external flow with freestream speed $V_{\infty}$:

$$ C_p \equiv \frac{p - p_{\infty}}{\frac{1}{2}\rho V_{\infty}^2} = 1 - \left(\frac{V}{V_{\infty}}\right)^2 $$