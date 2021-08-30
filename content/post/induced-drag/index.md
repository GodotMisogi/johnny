---
date: "2016-12-26T01:55:24+05:30"
title: "Calculus of Variations - Induced Drag Over a Wing"
draft: false
slug: induced-drag
description: A desperate attempt to come up with something original in aerodynamics using studies from physics, only having found it to be already discovered in the 1960s.
image: AircraftStream.svg
categories: 
    - Aerodynamics
    - Mathematics
---

While reading through John D. Anderson Jr.'s derivation of minimum induced drag, I thought of a cool application of the calculus of variations in one of the equations to deduce the required condition.

The equation that determines the downwash at a point is:

$$w(y\_0) = -\frac{1}{4\pi }\int^{b/2}\_{-b/2} \frac{(\mathrm{d}\Gamma/\mathrm{d}y)}{y\_0 - y}\mathrm{d}y = \int^{b/2}\_{-b/2} \mathcal{L}(\Gamma,\Gamma',y)~\mathrm{d}y$$ 

This effectively implies that the downwash can be expressed as a *functional* of $\Gamma$, i.e. $w\left[\Gamma(y)\right]$, and one can find the functional derivative to find the extremal point. There also exists a constraint on this system, the total lift across the span must be constant:

$$ L = \rho\_{\infty} V\_{\infty}\int^{b/2}\_{-b/2} \Gamma(y)~\mathrm{d}y = \int^{b/2}\_{-b/2} \mathcal{G}(\Gamma,\Gamma',y)~\mathrm{d}y$$

The Euler-Lagrange equations thus take the following form:

$$ \frac{\partial{\mathcal{L}}}{\partial{\Gamma}} - \frac{\mathrm{d}}{\mathrm{d}y}\left(\frac{\partial{\mathcal{L}}}{\partial{\Gamma'}}\right) + \lambda\left[\frac{\partial{\mathcal{G}}}{\partial{\Gamma}} - \frac{\mathrm{d}}{\mathrm{d}y}\left(\frac{\partial{\mathcal{G}}}{\partial{\Gamma'}}\right)\right]= 0 $$

Substituting the expressions:

$$ -\frac{1}{4\pi  (y\_0 - y)^2} + \rho\_{\infty}V\_{\infty}\lambda = 0$$

This doesn't contain any useful information about the downwash. Let's try something else.

Trying to minimise the induced drag formula directly as given by Anderson:

$$ C\_{D,i} = \frac{2}{V\_{\infty}S}\int^{b/2}\_{-b/2} \Gamma(x)\alpha\_{i}(x)~\mathrm{d}x = \frac{1}{2\pi V\_{\infty}^2 S}\int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\frac{\Gamma(x)\Gamma'(y)}{x - y}~\mathrm{d}y~\mathrm{d}x $$

Getting rid of the constants and performing a variation on the coefficient of induced drag, we get:

$$ \delta C\_{D,i} = \int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\left(\delta\Gamma(x)\frac{\Gamma'(y)}{x - y} + \delta\Gamma'(y)\frac{\Gamma(x)}{x - y}\right) ~\mathrm{d}y~\mathrm{d}x $$

Performing integration by parts on the second expression, keeping in mind that the first term of the evaluation disappears because the circulation at the endpoints (the boundary conditions of this problem) is zero:

$$ = \int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\delta\Gamma(x)\frac{\Gamma'(y)}{x - y}~\mathrm{d}y~\mathrm{d}x - \int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\delta\Gamma(y)\cdot\frac{\mathrm{d}}{\mathrm{d}y}\left(\frac{\Gamma(x)}{x-y}\right)~\mathrm{d}y~\mathrm{d}x$$

A little rearranging provides the more useful form:

$$ = \int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\delta\Gamma(x)\frac{\Gamma'(y)}{x - y}~\mathrm{d}y~\mathrm{d}x + \int^{b/2}\_{-b/2}\delta\Gamma(y)\cdot\frac{\mathrm{d}}{\mathrm{d}y}\int^{b/2}\_{-b/2}\frac{\Gamma(x)}{y-x}~\mathrm{d}x~\mathrm{d}y$$

A change of variables $y-x = q$ is useful to evaluate the last integral:

$$ \frac{\mathrm{d}}{\mathrm{d}y}\int^{b/2}\_{-b/2}\frac{\Gamma(x)}{y-x}~\mathrm{d}x = \frac{\mathrm{d}}{\mathrm{d}y}\int^{y-b/2}\_{y+b/2}\frac{\Gamma(y-q)}{q}~\mathrm{d}q$$

Feynman's favourite trick, differentiating under the integral sign:

$$ \require{cancel} \frac{\mathrm{d}}{\mathrm{d}y}\int^{y-b/2}\_{y+b/2}\frac{\Gamma(y-q)}{q}~\mathrm{d}q = \cancel{\frac{\Gamma(b/2)}{y-b/2}} - \cancel{\frac{\Gamma(-b/2)}{y+b/2}} + \int^{y-b/2}\_{y+b/2}\frac{\partial}{\partial y}\frac{\Gamma(y-q)}{q}~\mathrm{d}q $$

Mapping the variables back:

$$ \int^{y-b/2}\_{y+b/2}\frac{\partial}{\partial y}\frac{\Gamma(y-q)}{q}~\mathrm{d}q = \int^{b/2}\_{-b/2}\frac{\Gamma'(x)}{y-x}~\mathrm{d}x$$

Substituting this into the original expression:

$$ = \int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\delta\Gamma(x)\frac{\Gamma'(y)}{x - y}~\mathrm{d}y~\mathrm{d}x + \int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\delta\Gamma(y)\frac{\Gamma'(x)}{y - x}~\mathrm{d}x~\mathrm{d}y $$

Switching the variables of integration in the second expression, we get:

$$\delta C\_{D,i} = 2\int^{b/2}\_{-b/2}\int^{b/2}\_{-b/2}\delta\Gamma(x)\frac{\Gamma'(y)}{x - y}~\mathrm{d}y~\mathrm{d}x $$

Reintroducing the constants and combining this with the constraint, $\delta C\_{D,i} - \lambda\delta L = 0 $ becomes:

$$ \int^{b/2}\_{-b/2}\delta\Gamma(x)~\mathrm{d}x \left[\int^{b/2}\_{-b/2}\frac{2\Gamma'(y)}{x - y}~\mathrm{d}y\ - {2\pi\lambda}\right] = 0 $$

Using the constraint on the lift across the wing, this results in:

$$ \int^{b/2}\_{-b/2}\frac{\Gamma'(y)}{x - y}~\mathrm{d}y\ = \pi\lambda $$

The first term is part of the integral from the downwash expression at the beginning of the post, indicating that the downwash across the lifting line for minimum induced drag is constant:

$$ w = -\frac{\lambda}{4} = w\_0 $$

The same result as seen in Anderson, more rigorously!