---
date: "2020-08-30T22:55:27+08:00"
title: Gradient-Based Optimisation — Discrete Direct and Adjoint Methods (In Progress)
draft: false
slug: opt-problems
description: Demystifying the hype of the Newton method for optimisations, with relevant implementations.
image: XDSM.svg
categories: 
    - Physics
    - Mathematics
    - Programming
    - Computational Fluid Dynamics
---

Optimisation problems are becoming a big meme, so let's address some of the general ideas and how to compute them for discretised systems of equations.

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

$$ \frac{\partial R^{ijk}}{\partial p^{lmn}} = \ldots$$

Let us naïvely solve the linear system, by using Julia's `\` function from the `LinearAlgebra` package.

```julia
function newton_solver!(R, p, num_iters = 3, α = 1.0)
    # Array to store errors
    ε = zeros(num_iters);

    # Newton iteration loop
    for i in 1:num_iters
        # Compute residuals
        R = compute_residuals!(R, p)

        # Compute Jacobian
        ∂R∂p = compute_grad_fwd(R, p)

        # Compute Newton step
        Δp = ∂R∂p \ -R

        # Update state with relaxation factor
        p .= newton_update!(p, Δp, α)

        # Error processing
        ε[i] = maximum(abs.(Δp))
        println("Newton step error: $(ε[i])")
    end
    R, p, ε
end
```

There are various numerical methods dedicated to solving linear systems, depending on the structure of the matrices. Special methods depending on the case can be applied, for example using Jacobian-vector products instead of constructing the matrices directly. In the case of computational fluid dynamics with structured grids, the matrix is usually regular, banded, and sparse. This is due to the adjacencies between elements being specified up to a given level of "depth", which is called a stencil. So certain algorithms speed up the computations by taking advantage of this structure.


## Optimisations

The objective function for a gradient-based optimisation problem with an equality constraint (read holonomic, `bollu`) $\mathbf R(\mathbf x, \mathbf p(\mathbf x)) = \mathbf 0$ requires evaluation of its total derivative. The total derivative can be expanded it in two ways via the chain rule:

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} + \frac{\partial \mathbf f}{\partial \mathbf R}\frac{\partial\mathbf R}{\partial\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} + \frac{\partial \mathbf f}{\partial \mathbf p}\frac{d\mathbf p}{d\mathbf x} $$

The equality constraint trivially satisfies the following identity:

$$ \frac{d\mathbf R}{d\mathbf x} = \frac{\partial \mathbf R}{\partial \mathbf x} + \frac{\partial \mathbf R}{\partial \mathbf p}\frac{d\mathbf p}{d\mathbf x} = 0 $$

So you can substitute this into the second form of the total derivative expression:

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} + \frac{\partial \mathbf f}{\partial \mathbf p}\left[-\left(\frac{\partial \mathbf R}{\partial\mathbf p}\right)^{-1} \frac{\partial \mathbf R}{\partial \mathbf x}\right] $$

So there are two ways to avoid the problem of computing the inverse of $\partial \mathbf R/\partial \mathbf p$, which are called the _direct_ and _adjoint_ methods.

### Direct Method

Let

$$ \boldsymbol\psi = \frac{d \mathbf p}{d\mathbf x} = -\left(\frac{\partial \mathbf R}{\partial\mathbf p}\right)^{-1} \frac{\partial \mathbf R}{\partial \mathbf x} $$

The direct method hence solves the following linear system by left multiplication of $\partial \mathbf R / \partial \mathbf p$:
$$ \frac{\partial \mathbf R}{\partial \mathbf p}\boldsymbol \psi = -\frac{\partial \mathbf R}{\partial \mathbf x}$$

Notice that this expression is almost identical to the Newton step, extended to total derivatives with respect to the design variables. Essentially, the right-hand-side (RHS) is a matrix whose number of columns increases linearly with the number of design variables. This corresponds to numerous linear equations, one for each column of the RHS corresponding to a design variable. The big Jacobian on the left-hand-side is already computed from the Newton step, and only the RHS needs to be computed. If you have the factorisation of the Jacobian from the Newton step, it could be very efficient to simply reuse the factorisation and back-substitute to solve the numerous equations if the RHS is not very expensive to compute.

```julia
function solve_direct(x, u, ∂R∂x, ∂R∂u)
    ∂R∂u_sq = reshape(∂R∂u, (length(u[:]), length(u[:])))
    reshape(hcat((∂R∂u_sq \ -(∂R∂x)[:,:,i][:] for i in eachindex(x))...), (size(u)..., length(x)))
end
```

Hence the total derivative in this form is expressed as:

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} + \frac{\partial \mathbf f}{\partial \mathbf p}\boldsymbol\psi $$

```julia
total_derivative_direct(∂f∂x, ψ, ∂f∂u) 
    = ∂f∂x + [ sum(∂f∂u * ψ[n]) for n in eachindex(∂f∂x) ]
```

### Adjoint Method

Let 

$$ \boldsymbol\phi^T = \left(\frac{\partial \mathbf f}{\partial \mathbf R}\right)^T = \frac{\partial \mathbf f}{\partial \mathbf p}\left(\frac{\partial \mathbf R}{\partial\mathbf p}\right)^{-1} $$

The adjoint method hence solves the following linear system by right multiplication of $\partial \mathbf R / \partial \mathbf p$:

$$\left(\frac{\partial \mathbf R}{\partial \mathbf p}\right)^T \boldsymbol\phi = \left(\frac{\partial \mathbf f}{\partial \mathbf p}\right)^T $$

The transposition operator $-^T\colon (\mathbb R^m \times \mathbb R^n) \to (\mathbb R^n \times \mathbb R^m)$ simply transposes the dimensions of a matrix or vector. Of course, here the linear system has to be square, so $m = n$.

Here is where the word _adjoint_ comes in. When you have a constraint expressed in the form of some nonlinear partial differential equation, the transpose operator is not applicable, and you have to actually apply the concept of adjoints. This concept exists at a much more abstract level in category theory and duality, which is the "next-best" thing when one is not able to find direct inverses. An example is the adjoint of the differential operator on an infinite-dimensional Hilbert space of functions.

In most cases, however, you would be dealing with some discretised version of the PDE, and the above will apply to the discretised equations, which are expressed in the form of a tensor. This setup is called the _discrete adjoint_ formulation, and is quite a general operation on nonlinear PDEs.

```julia
function solve_adjoint(u, ∂R∂u, dfdu) 
    reshape(∂R∂u, (length(u[:]), length(u[:])))' \ (dfdu)'[:]
end
```

Hence the total derivative in this form is expressed as:

$$ \frac{d\mathbf f}{d\mathbf x} = \frac{\partial \mathbf f}{\partial \mathbf x} - \boldsymbol\phi^T\frac{\partial \mathbf R}{\partial \mathbf x} $$

```julia
total_derivative_adjoint(∂f∂x, φ, ∂R∂x) = 
    ∂f∂x + [ sum(permutedims(φ) * reshape(∂R∂x, (length(R[:]), length(∂f∂x)))[:,n]) 
             for n in eachindex(∂f∂x) ]
```

### Concept

So the main idea is that depending on the numbers of inputs and outputs, you choose the appropriate method which reduces the dimensionality of the linear system. When the number of outputs is larger than the number of inputs, use the direct method, and conversely use the adjoint method. The complicated parts (read headaches) are actually the evaluations and constructions of the different Jacobians and their reshapings, all of which are needed regardless of whether the direct or adjoint method is used. You picks your linear system and you takes your computational cost.

Notice that we didn't specify any restriction on the constraints. You could identically apply the same setup to a system of nonlinear equations:

$$\mathbf R(\mathbf x, \mathbf p_1, \mathbf p_2, \ldots, \mathbf p_n) = \begin{bmatrix} \mathbf R_1(\mathbf x, \mathbf p_1) \\\\ \mathbf R_2(\mathbf x, \mathbf p_2) \\\\ \vdots \\\\ \mathbf R_n(\mathbf x, \mathbf p_n) \end{bmatrix} = \mathbf 0$$

This setup is called a _coupled_ system, and the entire system is solved simultaneously for robustness and efficiency, either using the direct or adjoint methods. The "off-diagonal" block terms correspond to connections between different residual equations, and usually make the system sparse as many variables are not connected between different residual equations in physical models. The complexity of such systems is much greater, and is better left for another post.

## Stability

Even the order of accuracy of the Jacobians is not as important, as the point of the optimisation is to get a converged optimum, with converged residuals for the constraint. So you could spend less time doing more iterations with a greatly reduced cost for evaluating less accurate Jacobians and solving the linear system, if the residual constraint is the most time-consuming part. However, robustness of the solver for the residual equations becomes problematic.

<!-- ### Pseudo-Transient Continuation

$$ \left(\frac{I}{\Delta t} + 1\right)? \mathbf J_i d\mathbf x_i = -\mathbf f_i $$ -->

<!-- ## Unsteady Problems -->