---
date: "2021-02-24T22:55:27+08:00"
title: Gradient-Based Optimization - Discrete Direct and Adjoint Methods
draft: true
slug: opt-problems
description: Demystifying the hype of the Newton method for optimisations, with relevant implementations.
image: XDSM.svg
categories: 
    - Physics
    - Mathematics
    - Programming
    - Computational Fluid Dynamics
---

Optimisation problems are becoming a big meme, so let's address some of the general ideas and how to compute them on discretised systems of equations.

## Problem Description

Let the vector $\mathbf p_i \in \mathbb R^n$ be the state variable at the $i$th iteration, and the function as a matrix $\mathbf R \in \mathbb R^m \times \mathbb R^n$ which needs to be solved.

$$ \mathbf R(\mathbf p) = \mathbf 0 $$

Now say you want to know the values of $\mathbf p$, i.e. a root, which solves the system $\mathbf R = \mathbf 0$. The first step would be to consider some point close to the root $\mathbf p + \Delta \mathbf p$, then the system can be expanded in a Taylor series.

$$ \mathbf R(\mathbf p) \approx \mathbf R(\Delta \mathbf p) + \frac{\partial \mathbf R}{\partial \mathbf p}(\Delta \mathbf p) + \mathcal O(\Delta \mathbf p^2)$$

## Newton Method

The Newton method with a line-search performs the following algorithm:

$$ \mathbf p_{i+1} = \mathbf p_i - \alpha_i\mathbf J_i^{-1} \mathbf R(\mathbf p_i), \quad \mathbf J_i = \frac{\partial \mathbf R(\mathbf p_i)}{\partial \mathbf p_i} $$

assuming the existence of the inverse of the Jacobian $\mathbf J_i$, which is usually the case for physical problems (if not, re-evaluate your setup!). We would like $\mathbf p_{i+1} = 0$ in the fewest number of iterations or at the lowest speed. Hence $\mathbf p_i = \alpha_i\mathbf J_i^{-1} \mathbf R(\mathbf p_i)$. Evaluating the inverse of the Jacobian is impractical, but the required step $d\mathbf p_i$ can be determined by solving the following linear equation instead, which is simply multiplying both sides by $\mathbf J_i$.

$$ \beta_i \mathbf J_i d\mathbf p_i = -\mathbf R_i, \quad \text{where} \quad \beta_i = \frac{1}{\alpha_i} $$

The dimensionality of the problems could get confusing, so index notation with summation convention comes in handy:

$$ \frac{\partial R_{ijk}}{\partial p_{mn}} = $$

Let us naïvely solve the linear system, by using Julia's `\` function from the `LinearAlgebra` package.

```julia
# Number of iterations
num_iters = 5

# Newton step computation
function newton_step(u, R, ∂R∂u)
    i, j = size(R)
    l, k = size(u)
    jac = reshape(∂R∂u, (i*k, j*l))
    reshape(jac \ -R[:], size(u)) 
end

# Newton method algorithm
for i in 1:num_iters
    R, ∂R∂T = compute_residual_and_grad(T_truth, βs, q₀, dx, dy)
    ΔT = newton_step(T_truth, R, ∂R∂T)
    ε[i] = maximum(abs.(ΔT))
    # display(R)
    println("Newton step error: $(ε[i])")
    T_truth .+= ω * ΔT
end
```

There are various numerical methods dedicated to solving linear systems, depending on the structure of the matrices. In computational fluid dynamics with structured grids, the matrix is usually regular, banded, and sparse. This is due to the adjacencies between elements being specified up to a given level of "depth", which is called a stencil. So certain algorithms speed up the computations by taking advantage of this structure.


## Optimisations

The objective function for a gradient-based optimisation problem with an equality constraint (read holonomic, `bollu`) $\mathbf R(\mathbf x, \mathbf p) = \mathbf 0$ requires evaluation of its total derivative.

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} + \frac{\partial \mathbf f}{\partial \mathbf p}\frac{d\mathbf p}{d\mathbf x} $$

The equality constraint trivially satisfies the following identity:

$$ \frac{d\mathbf R}{d\mathbf x} = \frac{\partial \mathbf R}{\partial \mathbf x} + \frac{\partial \mathbf R}{\partial \mathbf p}\frac{d\mathbf p}{d\mathbf x} = 0 $$

So you can substitute this into the total derivative expression:

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} + \frac{\partial \mathbf f}{\partial \mathbf p}\left[-\left(\frac{\partial \mathbf R}{\partial\mathbf p}\right)^{-1} \frac{\partial \mathbf R}{\partial \mathbf x}\right] $$

So there are two ways to avoid the problem of computing the inverse of $\partial \mathbf R/\partial \mathbf P$, which are called the _direct_ and _adjoint_ methods.

### Direct Method

Let

$$ \boldsymbol\psi = -\left(\frac{\partial \mathbf R}{\partial\mathbf p}\right)^{-1} \frac{\partial \mathbf R}{\partial \mathbf x} $$

The direct method hence solves the following linear system by left multiplication of $\partial \mathbf R / \partial \mathbf p$:
$$
\begin{aligned}
    \frac{\partial \mathbf R}{\partial \mathbf p}\boldsymbol \psi = -\frac{\partial \mathbf R}{\partial \mathbf x}, \quad \text{where} \quad \boldsymbol \psi = \frac{d \mathbf p}{d\mathbf x}
\end{aligned}
$$

Notice that this expression is almost identical to the Newton step, extended to total derivatives with respect to the design variables. The big Jacobian on the LHS is already computed from the Newton step, and only the RHS needs to be computed. If you have the factorisation of the Jacobian from the Newton step, it could be very efficient to simply reuse the factorisation and back-substitute to solve this equation.

```julia
function solve_direct(x, u, ∂R∂x, ∂R∂u)
    ∂R∂u_sq = reshape(∂R∂u, (length(u[:]), length(u[:])))
    reshape(hcat((∂R∂u_sq \ -(∂R∂x)[:,:,i][:] for i in eachindex(x))...), (size(u)..., length(x)))
end
```

Hence the total derivative is:

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} + \frac{\partial \mathbf f}{\partial \mathbf p}\boldsymbol\psi $$

```julia
total_derivative_direct(∂f∂x, ψ, ∂f∂u) 
    = ∂f∂x + [ sum(∂f∂u * ψ[n]) for n in eachindex(∂f∂x) ]
```

### Adjoint Method

Let 

$$ \boldsymbol\phi^T = \frac{\partial \mathbf f}{\partial \mathbf p}\left(\frac{\partial \mathbf R}{\partial\mathbf p}\right)^{-1} $$

The adjoint method hence solves the following linear system by right multiplication of $\partial \mathbf R / \partial \mathbf p$:

$$
\begin{aligned}
    \left(\frac{\partial \mathbf R}{\partial \mathbf p}\right)^T \boldsymbol\phi = \left(\frac{\partial \mathbf f}{\partial \mathbf p}\right)^T, \quad \text{where}\quad \boldsymbol \phi = \left(\frac{\partial \mathbf f}{\partial \mathbf R}\right)^T
\end{aligned}
$$

The transposition operator $-^T\colon (\mathbb R^m \times \mathbb R^n) \to (\mathbb R^n \times \mathbb R^m)$ simply transposes the dimensions of a matrix or vector. Of course, here the linear system has to be square, so $m = n$.

Here is where the word _adjoint_ comes in. When you have a constraint expressed in the form of some nonlinear partial differential equation, the transpose operator is not applicable, and you have to actually apply the concept of adjoints. This concept exists at a much more abstract level in category theory and duality, in which one is not able to find direct inverses? An example is the adjoint of the differential operator on an infinite-dimensional Hilbert space of functions.

In most cases, however, you would be dealing with some discretised version of the PDE, and the above will apply to the discretised equations, which are expressed in the form of a matrix. This setup is called the _discrete adjoint_ formulation, and is quite a general operation on nonlinear PDEs.

```julia
function solve_adjoint(u, ∂R∂u, dfdu) 
    reshape(∂R∂u, (length(u[:]), length(u[:])))' \ -(dfdu)'[:]
end
```

Hence the total derivative is:

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} - \boldsymbol\phi^T\frac{\partial \mathbf R}{\partial \mathbf x} $$

```julia
total_derivative_adjoint(∂f∂x, φ, ∂R∂x) 
    = ∂f∂x + [ sum(permutedims(φ) * 
                reshape(∂R∂x, (length(R[:]), length(∂f∂x)))[:,n]) 
                for n in eachindex(∂f∂x) ]
```

### Concept

So the main idea is that depending on the number of inputs and outputs, you choose the appropriate method which reduces the dimensionality of the linear system. 

The complicated parts (read headaches) are actually the evaluations and constructions of the different Jacobians and their reshapings, all of which are needed regardless of whether the direct or adjoint method is used. You picks your linear system and you takes your computational cost.

Notice that we didn't specify any restriction on the constraints. You could identically apply the same setup to $$\mathbf R(\mathbf x, \mathbf p_1, \mathbf p_2, \ldots, \mathbf p_n) = \begin{bmatrix} \mathbf R_1(\mathbf x, \mathbf p_1) \\\\ \mathbf R_2(\mathbf x, \mathbf p_2) \\\\ \vdots \\\\ \mathbf R_n(\mathbf x, \mathbf p_n) \end{bmatrix} = \mathbf 0$$

This setup is called a _coupled_ system, and the entire system is solved simultaneously for robustness and efficiency, either using the direct or adjoint methods. The "off-diagonal" block terms correspond to connections between different residual equations, and usually make the system sparse as many variables are not connected between different residual equations in physical models. The complexity of such systems is much greater, and is better left for another post.

## Stability

Even the order of accuracy of the Jacobians is not as important, as the point of the optimisation is to get a converged optimum, with converged residuals for the constraint. So you could spend less time doing more iterations with a greatly reduced cost for evaluating less accurate Jacobians and solving the linear system, if the residual constraint is the most time-consuming part. However, robustness of the solver for the residual equations becomes problematic.

### Pseudo-Transient Continuation

$$ \left(\frac{I}{\Delta t} + 1\right)? \mathbf J_i d\mathbf x_i = -\mathbf f_i $$
