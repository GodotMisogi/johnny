---
date: "2021-09-03T12:54:58+08:00"
title: "Category Theory — Back to School"
draft: true
description: Category theory is a progressive scam for re-studying arithmetic, algebra, calculus, and category theory.
slug: category-theory
image: 
categories:
    - Mathematics
    - Programming
---

## Introduction

I studied category theory with `Vyn` from Bartosz Milewski's book/blog _Category Theory For Programmers_ in 2020, after a brief introduction from my course in _Advanced Algebra_ at HKUST. Eventually we realised it's just a scam to study mathematical concepts starting from primary school to itself. I list the concepts from memory which I usually find difficult to recall here, with attempted proofs without reference. The concepts are presented in the temporal sequence of standard mathematical education.

Some notation:

1. Categories will be notated using bold with a capital, e.g. $\textbf{Mon}$ for the category of monoids.
2. Objects and morphisms of a category $\mathbf C$ will be denoted as $\mathrm{Obj}(\mathbf C),~\mathrm{Mor}(\mathbf C)$ respectively.
3. Functor $\mathcal F$-mappings of objects $A \in \mathrm{Obj}(\mathbf C)$ will be denoted with curly brackets, $\mathcal F(A)$, and of morphisms $f \in \mathrm{Mor}(\mathbf C)$ with square brackets, $\mathcal F[f]$.
4. The category of functors between two categories will be notated in brackets, e.g. $[\mathbf{Grp},\mathbf{Set}]$ for the category of functors from the category of groups to the category of sets.

Milewski provides a nice programmatic interpretation of functors as, to paraphrase, "contextual wrappers" around morphisms between objects. Haskell examples are useful for illustrative purposes, and have been presented whenever appropriate to elucidate this idea.

## Arithmetic

### Yoneda Lemma

Every preschooler has learnt the Yoneda lemma in kindergarten. Here we recall the general concepts and proofs.

**Definition.** A locally small category is a category in which the morphisms between any two objects form a set.

Let $\mathbf C$ be a locally small category. Let $A \in \mathrm{Obj}(\mathbf C)$, $B \in \mathrm{Obj}(\mathbf{Set})$.

Denote the set-valued hom-functor as $\mathbf C(A, -)\colon \mathbf C \to \mathbf{Set}$ and let $\mathcal F\colon \mathbf C \to \mathbf{Set}$ be a covariant functor. 

A natural transformation between the two functors $\alpha \colon \mathbf C(A,-) \to \mathbf{Set}$ is expressed as a set of maps $\alpha_A, \alpha_B$ for all $A \in \mathrm{Obj}(\mathbf{C}), B \in \mathrm{Obj}(\mathbf{Set})$.

_Lemma_. (Yoneda) The set of natural transformations $\alpha$ is in one-to-one correspondence with $\mathcal F(A)$. Similarly for contravariant $\mathcal G\colon \mathbf C \to \mathbf{Set}$.

_Proof_. 

<!-- $$\begin{CD} \mathcal C(A, A) @>{\alpha_A}>> \mathcal C(A, B)\end{CD}$$ -->

```
Hom(A, A) —Hom(A, f)→ Hom(A, B)S
    |                     |
   α_A                   α_B
    ↓                     ↓
   C(A)   ----C(f)--→    C(B)
```


**Example.**

### Universal Properties

The formal definition of a universal property is usually provided in terms of a comma category, or shape categories and diagram functors (shown later). Here we look at some special definitions which are "arithmetic" in nature.

#### Products and Coproducts (Sums)

A product is defined by the following diagram:

**Example.** (Product) Cartesian product of two sets?

A product is defined by the following diagram, which applies a complicated technique called _dualisation_ to the product diagram, viz. flipping the arrows:

**Example.** (Coproduct) Disjoint union of two sets?

In Haskell, the coproduct is represented as a sum type with the "divider" syntax `|`, for example with the `Maybe` and `Either` types:

```haskell
data Maybe a = Just a | Nothing
data Either a b = Left a | Right b
```

> **Hence we have been scammed to re-learn arithmetic.**

## Algebra

### Monoidal Categories

In algebra, a monoid is defined by an associative multiplication $\mu$ and a unit map $\eta$. A monoidal category defines a unit map over a monoidal product (canonically denoted using the tensor product symbol $\otimes$), and assumes the existence of a special object $I \in \mathbf C$ such that a unit map $\eta \colon I \to A$ is defined for any object $A \in \mathrm{Obj}(\mathbf C)$. The multiplication, unit, and associators are encoded in the following commutative diagrams:

### Monads

In the words of Saunders Mac Lane: "A monad is just a monoid in the category of endofunctors". This (intentionally obfuscating) definition is ubiquitous for scaring programmers, and perhaps anyone unfamiliar with category theory. Here we attack the problem directly using the definitions.

**Example.** (List monad)

`List` is the canonical example of a monad as a polymorphic recursive sum type.
```haskell
data List a = Nil | Cons a (List a)
```

First we define the functor-mapping rules:
```haskell
instance Functor List where
    fmap f Nil         = Nil
    fmap f (Cons y ys) = Cons (f y) (fmap f ys) 
```

And now we define the monoidal laws for appending two lists by recursion:
```haskell
instance Monoid List where
    mappend :: List a -> List a -> List a
    mappend ys Nil         = ys
    mappend Nil zs         = zs
    mappend (Cons y ys) zs = Cons y (mappend ys zs)
    
    mempty :: List a
    mempty = Nil
```

The monoidal multiplication $M \otimes M \to M \in \mathrm{Mor}(\mathbf C)$. Now consider the category of endofunctors $[\mathbf C, \mathbf C]$ as a monoidal category, in which the monoidal product is composition of the endofunctors $\mathcal T \circ \mathcal T \to \mathcal T \in \mathrm{Mor}([\mathbf C, \mathbf C])$.

The multiplication map is list concatenation, and the unit map is creation of a single-element list from the argument. To define the monad, we just need to use the monoidal laws with composition and recursion.

```haskell
instance Monad List where
    join :: List (List a) -> List a
    join Nil           = mempty
    join (Cons ys zss) = mappend ys (join zss)

    return :: a -> List a
    return y = Cons y mempty
```

It is conceptually simpler to understand the use of a monad in terms of a Kleisli category $\mathbf C^{\mathcal T}$ which is constructed using an endofunctor $\mathcal T \in [\mathbf C, \mathbf C]$, where the objects are $\mathrm{Obj}(\mathbf C^\mathcal T) = \mathrm{Obj}(\mathbf C)$ and morphisms are $A \to \mathcal T(B) \in \mathrm{Mor}(\mathbf C^\mathcal T)$ for $A,B \in \mathrm{Obj}(\mathbf C)$.

Monads are extremely useful for separating pure computation from impure work (such as logging or any other operations involving side-effects). Following this style in any language is very useful for debugging and writing error-free code.

**Example.** (State monad?)

```haskell
data State s a = State (s -> a) s
```

### Equivalence of Adjunctions and Unit-Counit Pairs

Two functors $\mathcal L\colon \mathbf C \to \mathbf D, \mathcal R \colon \mathbf D \to \mathbf C$ are called _adjoint_ if they satisfy the following relations on hom-functors:

$$\mathbf D(\mathcal L~A, B) \cong \mathbf C(A, \mathcal R~B)$$

Specifically, $\mathcal L$ is _left-adjoint_ to $\mathcal R$ and $\mathcal R$ is _right-adjoint_ to $\mathcal L$. The pair $(\mathcal L, \mathcal R)$ is called an _adjunction_.

A _unit_ is defined as a natural transformation $\eta \colon \Delta_I \to \mathcal{I}$ with $\mathcal L \circ \mathcal R \circ \mathcal L \to \mathcal I$ and a _counit_ is a natural transformation $\varepsilon\colon \mathcal I \to \mathcal R \circ \mathcal L \circ \mathcal R$.

**Example.** (Free-forgetful adjunction?)

### Algebras and Coalgebras

Algebras and coalgebras are special categories defined using monads, with definitions as maps $\mathcal T(a) \to a$ and $a \to \mathcal T(a)$ respectfully. Certain categories are defined by using the following diagrams for objects and morphisms:

**Example.** (Fold catamorphism for list monad)

```haskell
foldl :: (a -> a -> a) -> a -> List a -> a
foldl f x Nil = x
foldl f x (Cons y ys) = foldl f (f x y) ys
```

Using this, the definition of the monadic law `join` for `List` reduces to a special case of applying the monoidal laws to `foldl`:
```haskell
join :: List (List a) -> List a
join = foldl mappend mempty
```

> **Hence we have been scammed to re-learn algebra.**

## Calculus

### Limits and Colimits

Let $\mathbf J$ be a category (with some arbitrary construction, usually referred to as a shape category), and $\mathbf C$ be a category (to be focused on). Then one defines a _diagram functor_ $\mathcal D\colon \mathbf J \to \mathbf C$. This allows the construction of universal properties such as (co)products.

### Profunctors

Profunctors are special cases of bifunctors, which are set-valued functors defined on a product category $\mathbf D \times \mathbf C$, such that they are necessarily contravariant on the first argument (by convention). Hence a profunctor $\mathcal F \colon \mathbf D^{\mathrm{op}} \times \mathbf C \to \mathbf{Set}$. Wedges describe morphism mappings of profunctors?

**Example.** (Hom-profunctor)

### Ends and Coends

The end is a terminal wedge, and the coend is an initial cowedge? These are denoted using integrals:

$$\begin{aligned} \mathrm{End_\mathbf C} & = \int_{A \in \mathrm{Obj}(\mathbf C)} \mathcal F (A, A) \\\\ \mathrm{Coend_\mathbf C} & = \int^{A \in \mathrm{Obj}(\mathbf C)} \mathcal F (A, A)
\end{aligned}$$

These may look complicated, but they are basically the equivalents of integrals, viz. the limit of a product or a sum.

(Surprisingly, a very interesting concept called a _Kan extension_ can be defined in terms of an end.)

> **Hence we have been scammed to re-learn calculus.**
<!-- 
## Number Theory

### Lawvere Theories

> **Hence we have been scammed to re-learn number theory.** -->

## Category Theory

### 2-Categories

A $2$-category is a generalisation of a category, which is defined by the following axioms:

1. Existence of $0$-morphisms (analogous to objects of a category).
2. Existence of $1$-morphisms between $0$-morphisms (analogous to morphisms of a category).
3. Existence of $2$-morphisms between $1$-morphisms (no analogous construction in a category).

> **Hence we have been scammed to re-learn category theory.**