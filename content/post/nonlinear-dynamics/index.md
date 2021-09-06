---
date: "2021-04-06T21:53:06+05:30"
title: Nonlinear Dynamics — A Game of Computation (In Progress)
draft: false
description: Playing around with nonlinear dynamics and populations in Julia.
slug: nonlinear-dynamics
categories: 
    - Mathematics
    - Programming
---

## Volterra Equation

### Introduction

The Volterra equation governs coupled population dynamics between two species under a predator-prey relationship in a given environment.

Consider $x$ as the number of rabbits, and $y$ as the number of foxes. This is expressed in the form of the following autonomous system of ordinary differential equations. 

$$ 
\begin{aligned} 
    \dot x & = \alpha x - \eta yx, & \alpha, \eta > 0 \\\\ 
    \dot y & = -\beta y + \lambda xy, & \beta, \lambda > 0 \\\\ 
\end{aligned}
$$

Physically, these equations can be read off in the following manner: The change in the population of rabbits depends on a rate $\alpha$ proportional to the current population, and decreases at a rate $\eta$ proportional to the populations of foxes and rabbits, and similarly for the change in the population of foxes.

### Critical Points

The critical points of the system are the points at which the right-hand-sides are simultaneously zero. In this case, it's simple to evaluate by hand.

$$ 
\begin{aligned} 
    x(\alpha - \eta y) & = 0 \\\\ 
    y(\lambda x - \beta) & = 0
\end{aligned} \\\\ 
\implies (x_c, y_c) = \left\\{ (0, 0), \left(\frac{\beta}{\lambda}, \frac{\alpha}{\eta}\right) \right\\}
$$

#### Stability

At $(x_c, y_c) = (0, 0)$, any perturbation in $x$ while keeping $y$ fixed results in an increase in $\dot x$, hence it is unstable in this direction, and vice versa for $y$ and $\dot y$. This results in the following phase? plot.

At $(x_c, y_c) = \left(\dfrac{\beta}{\lambda}, \dfrac{\alpha}{\eta}\right)$, we must evaluate the system by linearisation.

### Computation

Let's do it ~~to~~ in Julia.
```julia
f1(x, y, α, η) =  α * x - η * y * x
f2(x, y, β, λ) = -β * y + λ * x * y
```

The `DifferentialEquations` package in Julia makes this very easy to evaluate. As this is an initial value problem, the following functions show the required setup:

```julia
function volterra!(R, xs, ps)
    R[1] = f1(xs[1], xs[2], ps[1,1], ps[1,2])
    R[2] = f2(xs[1], xs[2], ps[2,1], ps[2,2])

    R
end

α, η, β, λ = 1, 2, 3, 4
ps = [ α η; 
       β λ ]
x0 = [ β/λ, α/η ]
```

We now set up a problem in the format of `DifferentialEquations` and evaluate it:

```julia
using DifferentialEquations

tspan = (0.0, 500.0)
prob = ODEProblem(volterra!, x0, tspan, ps)
sol = solve(prob, maxiters = 500)
```

Very short code! 

(_My second personality_: But probably not shorter than the Haskell equivalent.)

(_My third personality to the second_: Weird flex, but ok.)

<center>
{{< load-plotly >}}
{{< plotly json="/post/nonlinear-dynamics/images/wabbitsfoxes.json" >}}
</center>

Notice how, after the maximum population of rabbits in a period, the population of foxes increases while the population of rabbits decreases to a local minimum. The foxes then quickly die off to a minimum and the population of rabbits increases, repeating the cycle. 

## Falkner-Skan Boundary Layer Equations

The Falkner-Skan transformation of the thin shear layer equations governing boundary layers is given by the following third-order differential equation:

$$ f''' + \left(\frac{1 + a}{2}\right) f f'' + a\left[1 - (f')^2\right] = 0 $$

with boundary conditions $f(0) = f'(0) = 0, ~ f'(\eta) = 1$, where $\eta$ is some specified upper bound. This can be expressed in terms of an autonomous system of first-order differential equations: $\dot{\mathbf x} = \mathbf f(\mathbf x, \mathbf p, η)$ by setting each derivative as an independent variable. So we define the function for the governing ODE: 

```julia
function falkner_skan_ODE!(du, x, p, η)
    F, U, S = x
    a       = p[1]
    du[1]   = U
    du[2]   = S             
    du[3]   = -(1 + a) / 2 * F * S - a * (1 - U^2)

    nothing
end
```

As this is a boundary value problem, in contrast to the previous initial value problem, another function is required to specify the boundary conditions. Note that the distinction between initial value and boundary value problems is "superficial", and the mathematical treatment is generally identical as the formulation of 'space' curves passing through specified points.

```julia
function falkner_skan_BC!(R, x, p, η)
    R[1] = x[1][1]              # η = 0, F = 0
    R[2] = x[1][2]              # η = 0, F' = 0
    R[3] = x[end][2] - 1        # η = ηₑ, F' = 1

    nothing
end
```

Now we feed the functions to `DifferentialEquations`, using a Mono-Implicit Runge-Kutta integration scheme to solve the discretised equations as a fully implicit system.

```julia
ηs    = (0.0, 10.0)
x0    = [5.0, 2.0, 0.0] 
as    = range(-0.2, 0.2, length = 10)

prob  = BVProblem.(Ref(falkner_skan_ODE!), Ref(falkner_skan_BC!), Ref(x0), Ref(ηs), as, syms = Ref([:y, :δ]))
sol   = solve.(prob, Ref(GeneralMIRK4()), dt = 2e-1)
```
<center>
{{< load-plotly >}}
{{< plotly json="/post/nonlinear-dynamics/images/lol.json" >}}
</center>

Here we can observe an inflection point at some specific value of $U/U_e$. Can we figure this out analytically?
