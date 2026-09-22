---
date: "2023-01-01T00:00:00+08:00"
title: "Abstract Algebra — Descending Endomorphisms"
draft: false
index: true
weight: 6
description: New constructs for mappings between group and graph theory.
image: descending-graph.svg
categories:
    - Research
    - Mathematics
---

**New constructs for mappings between group and graph theory.**

* Investigated a new definition for certain types of endomorphisms of groups called "descending" endomorphisms.
* Implemented programs in SageMath to compute descending endomorphisms of groups for proof validation.
* Analyzed representations of descending endomorphisms in graph theory.

A *descending* endomorphism $\delta$ of a group $G$ is one that descends to every quotient: for each normal subgroup $N \trianglelefteq G$, the quotient map $\varphi : G \to G/N$ satisfies $\varphi \circ \delta = \bar{\delta} \circ \varphi$ for a well-defined induced map $\bar{\delta}$ on $G/N$. These maps can be represented as directed graphs on the group's elements, connecting group and graph theory.

![Descending endomorphism graph of the dicyclic group Dic₄](descending-graph.svg)

### Publications

* Vinay Madhusudanan, Arjit Seth, and G. Sudhakara. "Descending Endomorphisms of Groups". *Palestine Journal of Mathematics* 12.1 (2023), pp. 318–325.
* Vinay Madhusudanan, Arjit Seth, and G. Sudhakara. "Descending endomorphisms of some families of groups". In: *Applied Linear Algebra, Probability and Statistics: A Volume in Honour of C. R. Rao and Arbind K. Lal*. Springer, 2023, pp. 409–424.
* Vinay Madhusudanan, G. Sudhakara, and Arjit Seth. "Descending endomorphism graphs of groups". *AKCE International Journal of Graphs and Combinatorics* 20.2 (2023), pp. 148–155. [doi:10.1080/09728600.2023.2234956](https://doi.org/10.1080/09728600.2023.2234956)
