---
date: "2021-04-06T21:53:06+05:30"
title: Nonlinear Dynamics — A Game of Computation (In Progress)
draft: false
description: Playing around with nonlinear dynamics and populations in Julia.
image: images/wabbitsfoxes.svg
slug: nonlinear-dynamics
categories: 
    - Mathematics
    - Programming
---

*Note:* The code presented over here is available on my repository [`NonlinearDynamics`](https://github.com/GodotMisogi/NonlinearDynamics), with the additional code that was used to generate the plots shown here.

## Lotka-Volterra Equation

### Introduction

The Lotka-Volterra equation governs coupled population dynamics between two species under a predator-prey relationship in a given environment. Here we will play with the populations of wabbits and foxes. Consider $x$ as the number of wabbits, and $y$ as the number of foxes. This is expressed in the form of the following autonomous system of ordinary differential equations, because populations are continuous variables.

$$ \begin{aligned} 
        \dot x & = \alpha x - \eta yx, & \alpha, \eta > 0 \\\\ 
        \dot y & = -\beta y + \lambda xy, & \beta, \lambda > 0 \\\\ 
    \end{aligned} $$

Physically, these equations can be read off in the following manner: The change in the population of wabbits depends on a rate $\alpha$ proportional to the current population, and decreases at a rate $\eta$ proportional to the populations of foxes and wabbits, and similarly for the change in the population of foxes.

### Critical Points

The critical points of the system $(x_c, y_c)$ are the points at which the right-hand-sides are simultaneously zero. In this case, it's simple to evaluate by hand.

$$ \begin{aligned} 
      x(\alpha - \eta y) & = 0 \\\\ 
      y(\lambda x - \beta) & = 0
    \end{aligned} \\\\ 
    \implies (x_c, y_c) = \left\\{ (0, 0), \left(\frac{\beta}{\lambda}, \frac{\alpha}{\eta}\right) \right\\} $$

#### Stability

At $(x_c, y_c) = (0, 0)$, any perturbation in $x$ while keeping $y$ fixed results in an increase in $\dot x$, hence it is unstable in this direction, and vice versa for $y$ and $\dot y$.

At $(x_c, y_c) = \left(\dfrac{\beta}{\lambda}, \dfrac{\alpha}{\eta}\right)$, we must evaluate the system by linearisation.

### Computation

Let's do it ~~to~~ in Julia. The [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/) package in Julia makes the setup and integration of this system very easy. The syntax of the code using the `@ode_def` macro from [`ParametrizedFunctions.jl`](https://github.com/SciML/ParameterizedFunctions.jl) looks very similar to the mathematical presentation of the system of first-order ODEs.

```julia
using DifferentialEquations

## ODE definition
lotka_volterra = @ode_def LotkaVolterra begin
    dx =  α * x - η * x * y    # Wabbit growth/decay
    dy = -β * y + λ * x * y    # Fox growth/decay
end α β η λ
```

As this is an initial value problem, we specify the initial conditions, parameters, and the timespan for evaluation. Let's start with $5$ wabbits and $1$ fox. The growth of wabbits is set to $2$ wabbits per second and the decay of foxes is set to $1$ fox per second. The decay of wabbits due to the fox population and the growth of foxes due to the rabbit population are both set to $1$ animal per second. Let's evaluate the evolution of this system for $1$ minute.

```julia
## Inputs and parameters
α, β  = 2, 1     # Growth/decay parameters
η, λ  = 1, 1     # Coupling parameters
ps    = [ α η ;  
          β λ ]
x0    = [5, 1]
tspan = (0.0, 60.0)
```

Now we prepare an `ODEProblem` and call `solve` to begin the integration of the differential equations. 
```julia
## ODEProblem and solution
run   = ODEFunction(lotka_volterra, syms = [:Wabbits, :Foxes])
prob  = ODEProblem(run, x0, tspan, ps)
sol   = solve(prob)
```

This setup is sufficient to run the problem. Very short code!  (_My second personality_: But probably not shorter than the Haskell equivalent.) (_My third personality to the second_: Weird flex, but ok.) Let's plot the results.

<center>
{{< load-plotly >}}
{{< plotly json="/post/nonlinear-dynamics/images/wabbitsfoxes.json" >}}
</center>

Notice how, after the maximum population of wabbits in a period, the population of foxes increases while the population of wabbits decreases to a local minimum. The foxes then quickly die off to a minimum and the population of wabbits increases, repeating the periodic cycle. 

Now let's observe the behaviour near the second critical point. As this point is stable, it is expected that the population dynamics will not vary much until the numbers reach the other critical point, so let's run it for a longer time of about $3$ minutes.

```julia
## Near stable critical point
x0    = [β / λ + 0.01, α / η + 0.01]
tspan = (0., 180.)
```

<center>
{{< load-plotly >}}
{{< plotly json="/post/nonlinear-dynamics/images/wabbitsfoxes-critical.json" >}}
</center>

Here we observe a periodic solution with low amplitude which appears to be slowly diverging as expected.

## Falkner-Skan Boundary Layer Equations

The Falkner-Skan transformation of the thin shear layer equations governing boundary layers is given by the following third-order differential equation:

$$ F''' + \left(\frac{1 + a}{2}\right) F F'' + a\left[1 - (F')^2\right] = 0 $$

with boundary conditions $F(0) = F'(0) = 0, ~ F'(\eta) = 1$, where $\eta$ is some specified upper bound. This can be expressed in terms of an autonomous system of first-order differential equations: $\dot{\mathbf x} = \mathbf f(\mathbf x, \mathbf p, η)$ by setting each derivative as an independent variable. 

$$ \frac{d}{dt} 
    \begin{bmatrix} 
      F \\\\ 
      U \\\\ 
      S 
    \end{bmatrix} = 
    \begin{bmatrix} 
      U \\\\ 
      S \\\\ 
      -\left(\frac{1 + a}{2}\right) FS - a \left(1 - U^2\right)
    \end{bmatrix} $$


### Computation

So we define the function for the governing ODE: 

```julia
falkner_skan = @ode_def FalknerSkan begin
    dF  = U                                   
    dU  = S                                   
    dS  = -(1 + a) / 2 * F * S - a * (1 - U^2)
end a
```

As this is a boundary value problem, in contrast to the previous initial value problem, another function is required to specify the boundary conditions. Note that the distinction between initial value and boundary value problems is "superficial", and the mathematical treatment is generally identical as the formulation of 'space' curves passing through specified points.

```julia
function falkner_skan_BC!(bcs, x, p, η)
    bcs[1] = x[1][1]        # η = 0,  F  = 0
    bcs[2] = x[1][2]        # η = 0,  F' = 0
    bcs[3] = x[end][2] - 1  # η = ηₑ, F' = 1

    nothing
end
```

Now we feed the functions as a `BVProblem` with the ODE and the boundary conditions, using a Mono-Implicit Runge-Kutta integration scheme to solve the discretised equations as a fully implicit system for stable integration.

```julia
ηs    = (0.0, 12.0)
x0    = [5.0, 2.0, 0.0] 
as    = sort([ range(-0.11, 0.2, length = 10); -0.092 ])

prob  = BVProblem.(Ref(falkner_skan_ODE!), Ref(falkner_skan_BC!), Ref(x0), Ref(ηs), as, syms = Ref([:y, :δ]))
sol   = solve.(prob, Ref(GeneralMIRK4()), dt = 2e-1)
```
<center>
{{< load-plotly >}}
{{< plotly json="/post/nonlinear-dynamics/images/falknerskan.json" >}}
</center>

Here we can observe the incepient separation point at a specific value of $a \approx -0.0904$. Can this be figured out analytically?

## Double Pendulum

Refer to [my previous post on the double pendulum](../../post/dubby-pendy/) for a review of the theory and equations. Julia's various libraries and their compositions allow me to set up the differential equations, choose the integration solver, and set up plots with layouts _very_ easily. This is much simpler than writing different differential equation solver algorithms and plotting tools manually for experimentation.

The double pendulum is a famous example of a chaotic system, and is very sensitive to time-step sizes, initial conditions, and even the solver algorithm itself due to compounding numerical noise.

### Stability and Fixed Points

The Lagrangian of the system is:

$$ \mathcal{L} = T - V  \quad \mathrm{where} \quad \begin{aligned} T & = \frac{1}{2}m_1 l_1^2 \dot{\theta}_1^2 + \frac{1}{2}m_2\left[l_1^2 \dot{\theta}_1^2 + l_2^2 \dot{\theta}_2^2 + 2l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos(\theta_1 - \theta_2)\right] \\\\ V & = -(m_1 + m_2)gl_1\cos \theta_1 - m_2gl_2\cos\theta_2 \end{aligned} $$

### Computation

Here the expressions are complicated and would benefit from temporary variables in the function, which are not supported in the domain-specific-language of `ParametrizedFunctions.jl`, so a regular function is prepared instead of using the `@ode_def` macro.

```julia
## ODE definition
function dubby_pendy!(dx, x, ps, t)
    θ₁, θ₂, p₁, p₂     = x
    l₁, l₂, m₁, m₂, g  = ps

    C₀  = l₁ * l₂ * (m₁ + m₂ * sin(θ₁ - θ₂)^2)
    C₁  = (p₁ * p₂ * sin(θ₁ - θ₂)) / C₀
    C₂  = (m₂ * ((l₂ * p₁)^2 + (m₁ + m₂) * (l₁ * p₂)^2) 
          - 2 * l₁ * l₂ * m₂ * p₁ * p₂ * cos(θ₁ - θ₂)) * sin(2(θ₁ - θ₂)) / (2C₀^2)

    dx[1] = (l₂ * p₁ - l₁ * p₂ * cos(θ₁ - θ₂)) / (l₁ * C₀)
    dx[2] = (l₁ * (m₁ + m₂) * p₂ - l₂ * m₂ * p₁ * cos(θ₁ - θ₂)) / (l₂ * m₂ * C₀)
    dx[3] = -(m₁ + m₂) * g * l₁ * sin(θ₁) - C₁ + C₂
    dx[4] = -m₂ * g * l₂ * sin(θ₂) + C₁ - C₂

    nothing
end 
```

Let's set up the initial conditions and parameters and solve the problem. For the purposes of smooth plotting with a "nice" trajectory, I've picked an adaptive Runge-Kutta method of order 4 solver with a small timestep; these choices directly influence the solutions themselves due to the chaotic nature of the system.

```julia
## Inputs and parameters
m1, m2 = 10.0, 12.0 # Masses
l1, l2 = 4.0, 7.0   # Lengths
g      = 9.81       # Gravitational acceleration

ps    = [ l1, l2, m1, m2, g ]
x0    = [ π/2, π, 1., -3. ]
tspan = (0.0, 60.0)
tstep = 1000
dt    = maximum(tspan) / tstep

## ODEProblem and solution
run   = ODEFunction(dubby_pendy!, syms = [:θ₁, :θ₂, :p₁, :p₂])
prob  = ODEProblem(run, x0, tspan, ps)
sol   = solve(prob, RK4(), dt = dt)
df    = DataFrame(sol)
```

<center>
<video autoplay loop><source src="images/dubbypendy.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
</center>
