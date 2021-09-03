---
date: "2021-09-03T12:54:58+08:00"
title: "Category Theory — Back to School"
draft: true
description: Category theory is a progressive scam for re-studying arithmetic, algebra, calculus, number theory, and category theory.
slug: category-theory
image: 
categories:
    - Mathematics
    - Programming
---

## Introduction

I studied category theory with `Vyn` from Bartosz Milewski's book/blog _Category Theory For Programmers_ in 2020, after a brief introduction from my course in _Advanced Algebra_ at HKUST. I list the concepts from memory which I usually find difficult to recall here, with attempted proofs without reference.

## Arithmetic

### Yoneda Lemma

Let $\mathbf C$ be a locally small category (a category in which the morphisms between any two objects form a set), and $\mathbf{Set}$ the category of sets. Let $A \in \mathbf C$, $B \in \mathbf{Set}$. 

Denote the set-valued hom-functor as $\mathbf C(A, -)\colon \mathbf C \to \mathbf{Set}$ and let $\mathcal F\colon \mathbf C \to \mathbf{Set}$ be a covariant functor. 

A natural transformation between the two functors $\alpha \colon \mathbf C(A,-) \to \mathbf{Set}$ is expressed as a set of maps $\alpha_A, \alpha_B$ for all $A \in \mathrm{Obj}(\mathbf{C}), B \in \mathrm{Obj}(\mathbf{Set})$.

_Lemma_. (Yoneda) The set of natural transformations $\alpha$ is in one-to-one correspondence with $F(A)$. Similarly for contravariant $\mathcal G\colon \mathbf C \to \mathbf{Set}$.

_Proof_. 

$$\begin{CD} \mathcal C(A, A) @>{\alpha_A}>> \mathcal C(A, B)\end{CD}$$

```
Hom(A, A) —Hom(A, f)→ Hom(A, B)S
    |                     |
   α_A                   α_B
    ↓                     ↓
   C(A)   ----C(f)--→    C(B)
```

### Universal Properties


## Algebra

### Monoidal Categories

In algebra, a monoid is defined by an associative multiplication $\mu$ and a unit map $\eta$. A monoidal category defines a unit map over a monoidal product (canonically denoted using the tensor product symbol $\otimes$), and assumes the existence of a special object $I \in \mathbf C$ such that a unit map $\eta \colon I \to A$ is defined for any object $A \in \mathbf C$.

### Monads

A monad is a monoid in the category of endofunctors $[\mathbf C, \mathbf C]$. 

### Equivalence of Adjunctions and Unit-Counit Pairs

Two functors $\mathcal L\colon \mathbf C \to \mathbf D, \mathcal R \colon \mathbf D \to \mathbf C$ are called _adjoint_ if they satisfy the following relations on hom-functors:

$$\mathbf D(\mathcal L~A, B) \cong \mathbf C(A, \mathcal R~B)$$

Specifically, $\mathcal L$ is _left-adjoint_ to $\mathcal R$ and $\mathcal R$ is _right-adjoint_ to $\mathcal L$. The pair $(\mathcal L, \mathcal R)$ is called an _adjunction_.

A _unit_ is defined as a natural transformation $\eta \colon \Delta_I \to \mathcal{I}$ with $\mathcal L \circ \mathcal R \circ \mathcal L \to \mathcal I$ and a _co-unit_ is a natural transformation $\varepsilon\colon \mathcal I \to \mathcal R \circ \mathcal L \circ \mathcal R$.

### Algebras and Coalgebras

## Calculus

### Limits and Colimits

Let $\mathbf J$ be a category (with some arbitrary construction, usually referred to as a shape category), and $\mathbf C$ be a category (to be focused on). Then one defines a _diagram functor_ $\mathcal D\colon \mathbf J \to \mathbf C$. This allows the construction of universal properties such as (co)products.

### Profunctors

### Ends and Coends

## Number Theory

### Lawvere Theories

## Category Theory

### 2-Categories