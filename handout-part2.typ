#import "@preview/cetz:0.5.2"
#import "@preview/lilaq:0.6.0" as lq
#import "@preview/physica:0.9.8": *

#import calc: even, odd

#set math.equation(numbering: "(1)")
#set page(numbering: "1")
#set heading(numbering: "1.")
#set par(justify: true)

#align(right)[Oliver Portee, 29.06.2026]
#align(center)[
#title[The Convergence Rate of a Multigrid Method with Gauss-Seidel Relaxation for the Poisson Equation]

*Part 2*
]


#let definition-box(body) = box(
	stroke: 1pt,
	inset: 7pt,
	width: 100%,
	body
)

#let fine-color = lq.color.map.petroff10.at(0)
#let coarse-color = lq.color.map.petroff10.at(1)

#let Gbox(content, color) = box(fill: color.transparentize(70%), inset: (top: 7pt, bottom: 4pt, left: 3pt, right: 2pt), radius: 5pt, stroke: color, baseline: 15%, content)
#let G1(content) = Gbox(content, fine-color)
#let G2(content) = Gbox(content, coarse-color)

#let W = 6
#let H = 4

#let grid-point(x, y) = {
	let fill = if even(x + y) {
		coarse-color
	} else {
		fine-color
	}
	cetz.draw.circle((x, y), fill: fill, stroke: none, radius: 4pt)
}


#let draw-coarse-grid() = {
	import cetz.draw: line
	for x in range(0, W + 1, step: 2) {
		line((x, 0), (x, H))
	}
	for y in range(0, H + 1, step: 2) {
		line((0, y), (W, y))
	}
	for x in range(W) {
		for y in range(H) {
			if odd(x + y) {
				line((x + 1, y), (x, y + 1))
			} else {
				line((x, y), (x + 1, y + 1))
			}
		}
	}
	for x in range(0, W + 1, step: 2) {
		for y in range(0, H + 1, step: 2) {
			grid-point(x, y)
		}
	}
	for x in range(1, W + 1, step: 2) {
		for y in range(1, H + 1, step: 2) {
			grid-point(x, y)
		}
	}
}

#let draw-th-grid() = {
	import cetz.draw: line

	for x in range(W + 1) {
		line((x, 0), (x, H))
	}
	for y in range(H + 1) {
		line((0, y), (W, y))
	}
	for x in range(W) {
		for y in range(H) {
			if odd(x + y) {
				line((x + 1, y), (x, y + 1))
			} else {
				line((x, y), (x + 1, y + 1))
			}
		}
	}
	for x in range(W + 1) {
		for y in range(H + 1) {
			if odd(x + y) {
				grid-point(x, y)
			}
		}
	}
}

#let draw-grid(show-th: false) = {
	import cetz.draw: line

	if show-th {
		let color = navy
		for x in range(1, W) {
			for y in range(1, H) {
				if odd(x + y) {
					line((x, y), (x + 1, y), (x, y + 1), close: true, fill: gradient.linear(color, white, white, angle: -45deg), stroke: none)
					line((x, y), (x + 1, y), (x, y - 1), close: true, fill: gradient.linear(color, white, white, angle: 45deg), stroke: none)
					line((x, y), (x - 1, y), (x, y - 1), close: true, fill: gradient.linear(color, white, white, angle: 135deg), stroke: none)
					line((x, y), (x - 1, y), (x, y + 1), close: true, fill: gradient.linear(color, white, white, angle: 225deg), stroke: none)
				}
			}
		}
	}
	
	for x in range(W + 1) {
		line((x, 0), (x, H))
	}
	for y in range(H + 1) {
		line((0, y), (W, y))
	}
	for x in range(W) {
		for y in range(H) {
			if odd(x + y) {
				line((x + 1, y), (x, y + 1))
			} else {
				line((x, y), (x + 1, y + 1))
			}
		}
	}
	for x in range(W + 1) {
		for y in range(H + 1) {
			grid-point(x, y)
		}
	}
}

= Introduction

The domain is discretized on level $q$ using the following grid.

#figure(
	cetz.canvas({
		cetz.draw.scale(0.8)
		draw-grid()
	}),
	caption: [Grid on level $q$. Yellow grid points belong to the coarser grid on level $q - 1$.]
)

The finite element space corresponding to this grid is $S_h$ (one degree of freedom at each inner grid point, piecewise linear, zero on the boundary). This function space can be written as a sum of two subspaces,
$
S_h = S_H plus.o T_h.
$
#figure(
	cetz.canvas({
		import cetz.draw: content, translate, scale
		scale(0.7)
		draw-grid()
		content((7, 2), $=$)
		translate(x: 8)
		draw-coarse-grid()
		content((7, 2), $plus.o$)
		translate(x: 8)
		draw-th-grid()
	}),
	caption: [The space $S_h$ of the fine grid can be expressed as a sum of the space $S_H$ on the coarse grid (middle) and the space $T_h$ (right).]
)

$T_h$ is defined as follows:
$
T_h := {w in S_h | w(p_i) = 0 "if" p_i "yellow"}.
$
The basis functions of $T_h$ form pyramid-like patterns. As a result, all functions in $T_h$ are zero on the diagonal lines through the coarse grid points.

#figure(
	cetz.canvas({
		cetz.draw.scale(0.8)
		draw-grid(show-th: true)
	}),
	caption: [Each shaded are corresponds to one of the basis functions of $T_h$. The supports of the respective basis functions are disjoint (they do not overlap each other).]
)


#let mg-scheme(
	u-labels: false,
	smoother-labels: false,
) = {
	import cetz.draw: line, content, scale, floating

	scale(1.3)
	let r = 4
	for x in range(r) {
		if x > 0 {
			line((x + 0.3, 0), (x + 0.7, 0), mark: (end: ">"))
		}
		if calc.odd(x) {
			content((x + 1, 0), G1($G^#(r - x)$))
		} else {
			content((x + 1, 0), G2($G^#(r - x)$))
		}
	}
	line((r + 0.1, -0.3), (r + 0.5, -1.0), mark: (end: ">"))
	content((r + 1, -1.5), align(center)[coarse grid correction])
	line((r + 1.5, -1), (r + 1.9, -0.3), mark: (end: ">"))
	for x in range(r + 1) {
		if x > 0 {
			line((r + 2 + x - 0.70, 0), (r + 2 + x - 0.3, 0), mark: (end: ">"))
		}
		if calc.even(x) {
			content((r + 2 + x, 0), G1($G^#(x + 1)$))
		} else {
			content((r + 2 + x, 0), G2($G^#(x + 1)$))
		}
	}

	if u-labels {
		content((0.5, -0.5), $u^0$)
		content((1.5, -0.5), $u^1$)
		content((2.5, -0.5), $u^2$)
		content((3.5, -0.5), $u^3$)
		content((4.5, -0.5), $u^r$)
		content((5.4, -0.5), $u^(r + 1)$)
		content((10.5, -0.5), $u^(2 r + 2)$)
	}

	if smoother-labels {
		cetz.decorations.brace((0.7, 0.4), (4.3, 0.4), name: "presmoother")
		content("presmoother.content", [presmoother $tilde(G)$])
		cetz.decorations.brace((5.7, 0.4), (10.3, 0.4), name: "postsmoother")
		content("postsmoother.content", [postsmoother $tilde(G)^*$])
	}
}

The multigrid algorithm consists of the presmoother, the coarse grid correction and the postsmoother.

#figure(
	cetz.canvas({
		mg-scheme(
			u-labels: true,
			smoother-labels: true,
		)
	}),
	caption: [Schematic representation of the multigrid algorithm for $r = 4$.#footnote[In the first iteration, we prepend $G^(r + 1)$ to the presmoother.]]
)

The individual steps in the smoothers are the Gauss-Seidel half-steps.
$
G^k := cases(G_h^"I" &wide& k "odd", G_h^"II" &wide& k "even")
$

The Gauss-Seidel half-steps are defined as follows:
$
(G_h^"I" u)_i &:= cases(u_i &wide& p_i "yellow", 1/4 sum_(j in N(i)) u_j + b_i &wide& p_i "blue") \
(G_h^"II" u)_i &:= cases(1/4 sum_(j in N(i)) u_j + b_i &wide& p_i "yellow", u_i &wide& p_i "blue")
$
Here, $N(i)$ is the set of the indices of the four neighboring grid points of $p_i$. The operator $G_h^"I" dot G_h^"II"$ represents a full Gauss-Seidel step. In the presmoother, we have $r$ half-steps,
$
u^r = G^r G^(r - 1) ... G^1 u^0 = tilde(G) u^0.
$
The only exception is iteration 1 (i.e., the first V-cycle or W-cycle) where we prepend the half-step $G^(r + 1)$ to achieve symmetry between presmoother and postsmoother. The postsmoother contains $r + 1$ half-steps,
$
u^(2 r + 2) = G^1 G^2 ... G^r G^(r + 1) u^(r + 1) = tilde(G)^* u^(r + 1)
$ <eq:Gstar>

Since we only want to examine convergence, we assume the homogeneous Laplace equation,
$
-laplace u &= 0 &wide& "in" Omega \
u &= 0 &wide& "on" partial Omega.
$
The exact solution of this equation is the zero function. Thus, $b_i = 0$ and a certain approximation $u^k$ is the same as the error with regard to the exact solution. Consequently, $G_h^"I"$ and $G_h^"II"$ simply compute the average of the function values at the four neighbors for each blue or yellow grid point.

We also introduce the energy functional
$
J(u) := a(u, u) - 2 (f, u)_0
$
where $a(dot, dot)$ is the bilinear form derived from the weak form of the poisson equation,
$
a(u, v) := integral_Omega grad u dot grad v = integral_Omega u_xi v_xi + u_eta v_eta,
$
and $(dot, dot)_0$ is the inner product of $L^2$,
$
(f, u)_0 := integral_Omega f u.
$
It can be shown that the solution of the strong (or weak) form of the poisson equation is equivalent to finding the minimizer of $J$. Note that for the homogeneous problem, $J$ simplifies to
$
J(u) = a(u, u).
$

Finally, we define the norm
$
||u|| := sqrt(a(u, u)).
$
By substituting the definition of $a(u, u)$, we can show that the norm can be computed via the following sum over pairs of grid points (remembering that each triangle has the area $h^2 / 2$):

$
||u||^2 = a(u, u) = sum_(T in cal(T)(Omega)) integral_(T) (u_"east" - u_"west")^2 / h^2 + (u_"north" - u_"south")^2 / h^2 = sum_vec(i\, j, d(i, j) = h, delim: #none) (u_i - u_j)^2.
$ <eq:norm>
Similarly, we define a semi-norm on the coarse grid points,
$
|u|^2 = sum_vec(i\, j "yellow", d(i, j) = 2 h, delim: #none) 1/2 (u_i - u_j)^2.
$ <eq:seminorm>
In @eq:norm and @eq:seminorm we only include unique pairs irrespective of the order (i.e., if we include the pair $(i, j)$ into the sum, then we do not include $(j, i)$).

In part 1 we have already shown the following three statements:

$
|a(v, w)| <= sqrt(1/2 (1 - (|v|^2) / (||v||^2))) ||v|| ||w|| wide v in S_H, w in T_h,
$ <eq:lemma-3-1>

$
|u| <= sqrt(rho) ||u|| wide "with" wide rho := (2 lambda^2) / (lambda^2 + 1) "and" lambda := (|v^r|) / (||v^r||),
$ <eq:3-11>

$
"if" wide u = G_h^"I" u wide "then" wide ||G_h^"II" u|| <= |u|.
$ <eq:lemma-4-1>

Now, we are looking for an upper bound of the convergence rate,
$
||u^(2 r + 2)|| <= delta_q ||u^0|| wide "with" wide delta_q < 1.
$
To this end, we first find the convergence rate of the presmoother (@sec:presmoother), then we include the coarse grid correction (@sec:exact-cgc and @sec:non-exact-cgc), and finally look at a full multigrid iteration (@sec:full-mg).

= Presmoother <sec:presmoother>

From @eq:3-11 and @eq:lemma-4-1, we find that if $u = G_h^"I" u$, then
$
||G_h^"II" u|| <= |u| <= sqrt(rho) ||u||.
$ <eq:G2-convergence>
Unfortunately, this only gives us an upper bound for the convergence rate of $G_h^"II"$ and not for $G_h^"I"$. Therefore, we will prove the following two statements
$
"if" u &= G_h^"II" u: &wide (||G_h^"I" u||) / (||u||) &<= (||G_h^"II" G_h^"I" u||) / (||G_h^"I" u||) \
"if" u &= G_h^"I" u: &wide (||G_h^"II" u||) / (||u||) &<= (||G_h^"I" G_h^"II" u||) / (||G_h^"II" u||)
$ <eq:ineffective>

First, we examine the $G_h^"I"$ in more detail.

#figure(
cetz.canvas({
	import cetz.draw: circle, line, content, scale, floating

	let grid-point(coords) = circle(
		coords,
		radius: 2pt,
		stroke: none,
		fill: if even(coords.at(0) + coords.at(1)) {
			fine-color
		} else {
			coarse-color
		}
	)

	let p0 = ( 0,  0)
	let p1 = ( 0, -1)
	let p2 = ( 1,  0)
	let p3 = ( 0,  1)
	let p4 = (-1,  0)
	let p5 = ( 1, -1)
	let p6 = ( 1,  1)
	let p7 = (-1,  1)
	let p8 = (-1, -1)

	scale(1.2)

	line(p5, p6, p7, p8, close: true)
	line(p1, p2, p3, p4, close: true)
	line(p1, p3)
	line(p2, p4)

	grid-point(p0)
	grid-point(p1)
	grid-point(p2)
	grid-point(p3)
	grid-point(p4)
	grid-point(p5)
	grid-point(p6)
	grid-point(p7)
	grid-point(p8)

	content(p0, $0$, anchor: "north-west", padding: 0.1)
	content(p1, $1$, anchor: "north", padding: 0.2)
	content(p2, $2$, anchor: "west", padding: 0.2)
	content(p3, $3$, anchor: "south", padding: 0.2)
	content(p4, $4$, anchor: "east", padding: 0.2)
}),
caption: [Neighborhood of a blue point with grid point numbering.]
) <fig:blue-neighborhood>

In this neighborhood, $G_h^"I"$ essentially computes the average of the four neighbors,
$
u_0 <- 1/4 (u_1 + u_2 + u_3 + u_4).
$
This operation minimizes the energy function $J(v)$ over the set $u + T_h$, since
$
J(v) = a(v, v) + cancel((f, v)_0) = ||v||^2 = sum_vec(i\, j, d(i, j) = h, delim: #none) (v_i - v_j)^2.
$ <eq:J-sum>
The part of the sum relevant for the neighborhood in @fig:blue-neighborhood is
$
(v_1 - v_0)^2 + (v_2 - v_0)^2 + (v_3 - v_0)^2 + (v_4 - v_0)^2.
$
This term can be interpreted as a function of $v_0$ (as $G_h^"I"$ only changes the values at the blue points). The minimum of this function can be computed by evaluating its derivative by $v_0$.
$
2 (v_0 - v_1) + 2 (v_0 - v_2) + 2 (v_0 - v_3) + 2 (v_0 - v_4) = 8 v_0 - 2 sum_(j = 1)^4 v_j =^! 0 \
==> v_0 = 1/4 (v_1 + v_2 + v_3 + v_4)
$
This is exactly what $G_h^"I"$ computes at each blue grid point, minimizing the sum in @eq:J-sum.

Now we define $v := G_h^"I" u$. Since $v$ is the minimizer of $J(u + w) quad forall w in T_h$, we can state that
$
J(v + w) >= J(v) wide forall w in T_h,
$ <eq:orthogonal-1>
or alternatively,
$
lr(dv(, epsilon) J(v + epsilon w)|)_(epsilon = 0) = 0 wide forall w in T_h, epsilon in RR.
$
Using the bilinear form $a(dot, dot)$, we expand this term as follows:
$
J(v + epsilon w) &= a(v + epsilon w, v+ epsilon w) \
&= a(v, v) + 2 epsilon a(v, w) + epsilon^2 a(w, w).
$
Then, the derivative is
$
dv(, epsilon) J(v + epsilon w) = 2 a(v, w) + 2 epsilon a(w, w).
$
For $epsilon = 0$, we get
$
lr(dv(, epsilon) J(v + epsilon w)|)_(epsilon = 0) = 2 a(v, w) =^! 0 wide forall w in T_h.
$ <eq:orthogonal-2>
Since $a(dot, dot)$ has the role of the inner product, we find that $v$ is perpendicular to $T_h$. If $T_h^perp := {u in S_h | (u, w) = 0 quad forall w in T_h}$ is the orthogonal complement of $T_h$, then $v in T_h^perp$.

Let $P_(T_h)$ and $P_(T_h^perp)$ denote the orthogonal projections onto $T_h$ and $T_h^perp$, respectively. $G_h^"I"$ only changes degrees of freedom that are part of $T_h$. Thus,
$
P_(T_h^perp) u = P_(T_h^perp) v = v = G_h^"I" u.
$
It follows that $G_h^"I"$ is the orthogonal projection onto $T_h^perp$.

This idea can also be visualized using a two-dimensional $S_h$ and one-dimensional $T_h$ and $T_h^perp$.

#let right-angle(pos, start: 0deg) = {
	import cetz.draw: arc, circle
	arc(pos, anchor: "origin", start: start, stop: start + 90deg)
	let r = 0.57
	let dx = r * calc.cos(start + 45deg)
	let dy = r * calc.sin(start + 45deg)
	let label-pos = (pos.at(0) + dx, pos.at(1) + dy)
	circle(label-pos, stroke: none, fill: black, radius: 1pt)
}

#figure(
	cetz.canvas({
		import cetz.draw: line, content, scale, floating, circle
		line((-1, -1), (4, 4), stroke: 1pt)
		content((4, 4), $T_h$, anchor: "west", padding: 0.3)
		content((4.5, 1.5), $v = G_h^"I" u in T_h^perp$, anchor: "south-west", padding: 0.2)
		content((3, 0), $u in S_h$, anchor: "north", padding: 0.2)
		content((1.5, 1.5), $P_(T_h) u = w in T_h$, anchor: "south-east", padding: 0.2)
		right-angle((3, 3), start: 225deg)
		circle((0, 0), radius: 3pt, stroke: none, fill: black)
		content((0, 0), $va(0)$, anchor: "south-east", padding: 0.1)
		line((0, 0), (6, 0), mark: (end: ">"), stroke: 2pt)
		line((0, 0), (3, 3), mark: (end: ">"), stroke: 2pt)
		line((3, 3), (6, 0), mark: (end: ">"), stroke: 2pt)
	}),
	caption: [Visualization of the operator $G_h^"I" = P_(T_h^perp)$.]
) <fig:G1>

Since $T_h plus.o T_h^perp = S_h$, we can find a $w in T_h$ such that $u - w = v$. Since $v perp w space forall w in T_h$, it follows that $w = P_(T_h) u$.

$G_h^"I"$ can also be expressed in terms of $P_(T_h)$:
$
v = u - w = u - P_(T_h) u = (I - P_(T_h)) u = G_h^"I" u,
$
where $I$ is the identity operator.

#highlight[The projection theorem tells us that the split of $u$ into $u = w + v$ is unique.]

Analogously, it can be shown that $G_h^"II"$ is an orthogonal projection. In this case, instead of $T_h$ we need to use the space $hat(T)_h := {w in S_h | w(p_i) = 0 "if" p_i "blue"}$

A fundamental property of orthogonal projections is that they are self-adjoint.

#definition-box[
	An operator $Q^*$ is the adjoint of $Q$ if
	$
	(Q u, v) = (u, Q^* v) wide forall u, v.
	$
	It is self-adjoint if $Q^* = Q$, i.e.,
	$
	(Q u, v) = (u, Q v) wide forall u, v.
	$
]

This property of orthogonal projections can be shown as follows. Let $x = x_1 + x_2 in S_h$ and $y = y_1 + y_2 in S_h$ such that $x_1, y_1 in T_h$ and $x_2, y_2 in T_h^perp$. If $P$ is the orthogonal projection onto $T_h$, then
$
(P x, y) = (P x_1 + cancel(P x_2), y) = (x_1, y) = (x_1, y_1 + y_2) \
= (x_1, y_1) + cancel((x_1, y_2)) = (x_1, y_1) = (x_1, y_1) + cancel((x_2, y_1)) \
= (x_1 + x_2, y_1) = (x, y_1) = (x, P y_1 + cancel(P y_2)) = (x, P y).
$ <eq:orthogonal-projections-self-adjoint>
The inner products $(x_1, y_2)$ and $(x_2, y_1)$ are zero because one of the arguments is in $T_h$ and the other in $T_h^perp$ (so they are orthogonal to each other).

Another property of orthogonal projections is that they are idempotent:
$
P^2 u = P P u = P u wide forall u.
$

Now we can prove @eq:ineffective. Assuming $u = G_h^"II" u$, we get
$
||G_h^"I" u||^2 &= (G_h^"I" u, G_h^"I" u) &wide& "definition of the norm" \
&= (u, G_h^"I" G_h^"I" u) &wide& G_h^"I" "self-adjoint" \
&= (u, G_h^"I" u) &wide& G_h^"I" "idempotent" \
&= (G_h^"II" u, G_h^"I" u) &wide& "by assumption" \
&= (u, G_h^"II" G_h^"I" u) &wide& G_h^"II" "self-adjoint" \
&<= ||u|| dot ||G_h^"II" G_h^"I" u|| &wide& "Cauchy-Schwarz inequality".
$ <eq:ineffective-1>

#definition-box[
Cauchy-Schwarz inequality:
$
|(u, v)| <= ||u|| dot ||v|| wide forall u, v
$ <eq:cauchy-schwarz>
The case of equality happens when $u$ and $v$ are linearly dependent,
$
|(u, v)| = ||u|| dot ||v|| wide <==> wide u = alpha v wide alpha in RR.
$
]
Rearranging @eq:ineffective-1, we get
$
"if" u = G_h^"II": wide (||G_h^"I" u||) / (||u||) <= (||G_h^"II" G_h^"I" u||) / (||G_h^"I" u||).
$
Similarly, this can be shown the other way around (when we swap $G_h^"I"$ and $G_h^"II"$):
$
"if" u = G_h^"I": wide (||G_h^"II" u||) / (||u||) <= (||G_h^"I" G_h^"II" u||) / (||G_h^"II" u||).
$
These two equations can be applied to the presmoother, which is an alternating sequence of $G_h^"I"$ and $G_h^"II"$:
$
(||u^(k + 1)||) / (||u^k||) <= (||u^(k + 2)||) / (||u^(k + 1)||) wide k = 0, ..., r - 2.
$ <eq:ineffective-2>
This means that the Gauss-Seidel half-steps become more ineffective over time (or in case of equality, the norm reduction stays the same). Note that the last step of the presmoother is always $G_h^"I"$, i.e., $u^r = G_h^"I" u^(r - 1)$. If we were to append another Gauss-Seidel half-step, we would have $u^(r + 1) = G_h^"II" u^r$. We know from @eq:G2-convergence that this step would have a norm reduction of $sqrt(rho)$ or less (the lower the better). Due to @eq:ineffective-2, his hypothetical last step $u^(r + 1) = G_h^"II" u^r$ is the most "ineffective" step:
$
(||u^(r + 1)||) / (||u^r||) >= (||u^(k + 1)||) / (||u^k||) wide k = 1, ..., r.
$
Therefore the first $r$ half-steps also have a norm reduction of $sqrt(rho)$ or smaller. Consequently, after $r$ half-steps, we have the following upper bound of the norm reduction:
$
||u^r|| <= (sqrt(rho))^r ||u^0|| = rho^(r / 2) ||u^0||.
$ <eq:presmoother-convergence>

= Exact Coarse Grid Correction <sec:exact-cgc>

First we look at an exact coarse grid correction (i.e., the equation is solved without error on the coarser grid). In this step of the multigrid algorithm, we minimize the energy function $J(v)$ over $u^r + S_H$. In other words, we evaluate
$
u^(q - 1) := limits("argmin")_(v in S_H) space J(u^r - v).
$
Here, $u^(q - 1)$ is in $S_H$, the finite element space associated with the coarser grid on level $q - 1$. In the correction step, we compute
$
u' = u^r - u^(q - 1).
$

Note that the exact coarse grid correction is similar to the Gauss-Seidel half steps, only that we minimize $J$ over $u^r + S_H$ instead of $u + T_h$. Hence, using the same argument as in @eq:orthogonal-1 to @eq:orthogonal-2, we find that $u' perp S_H$ and
$
u^(q - 1) = P_(S_H) u^r \
u' = u^r - u^(q - 1) = I u^r - P_(S_H) u^r = (I - P_(S_H)) u^r = P_(S_H^perp) u^r.
$ <eq:exact-cgc>

Furthermore, as observed earlier already, the last step of the presmoother is always $G_h^"I"$. The resulting vector $u^r$ is thus orthogonal to $T_h$ (compare vector $v$ in @fig:G1).

#figure(
	cetz.canvas({
		import cetz.draw: *
		
		scale(0.8)
		line((0, 0), (-4, 5), mark: (end: ">"), stroke: 2pt)
		content((-2, 2.5), $u^r$, anchor: "north-east", padding: 0.1)
		line((0, 0), (0, 5), mark: (end: ">"), stroke: 2pt)
		content((0, 2.5), $u'$, anchor: "west", padding: 0.2)
		line((0, 0), (25/4, 5), mark: (end: ">"), stroke: 2pt)
		content((25/8, 2.5), $w$, anchor: "north-west", padding: 0.1)
		line((0, 5), (-4, 5), mark: (end: ">"), stroke: 2pt)
		content((-2, 5), $u^(q - 1)$, anchor: "south", padding: 0.1)

		line((-7, 0), (7, 0))
		content((3.4, 0), $S_H$, anchor: "north", padding: 0.2)
		line((-7, 5), (7, 5), stroke: (dash: "dashed"))

		line((0.8, -1), (0, 0), stroke: (dash: "dashed"))

		right-angle((0, 0), start: 219deg)

		line((-1, -0.8), (15/2, 6), stroke: (dash: "dashed"))
		content((15/2, 6), $T_h$, anchor: "south-east", padding: 0.1)

		arc((0, 0), anchor: "origin", start: 90deg, stop: 128deg)
		content((-0.21, 0.7), $alpha$)
		arc((0, 0), anchor: "origin", start: 0deg, stop: 38deg)
		content((0.7, 0.29), $alpha$)

		right-angle((0, 5), start: 180deg)
	}),
	caption: [Visualization of the exact coarse grid correction, which is an orthogonal projection of $u^r$ onto $S_H^perp$.]
)

From @eq:lemma-3-1, we find
$
((v, w)) / (||v|| ||w||) = cos(alpha) <= sqrt(1/2 (1 - lambda^2)) wide forall v in S_H, w in T_h.
$
Since this equation is true for all $v in S_H$ and all $w in T_h$, this really says something about the angle $alpha$ between subspaces $S_H$ and $T_h$. Furthermore, this angle is the same as the angle between $u'$ and $u^r$ (as a consequence of $u' perp S_H$ and $u^r perp T_h$).

Due to the right angle between $u'$ and $u^(q - 1)$, we find
$
||u'|| = cos(alpha) ||u^r|| &<= sqrt(1/2 (1 - lambda^2)) ||u^r|| = sqrt(((1 - lambda^2) / (lambda^2 + 1)) / (2 / (lambda^2 + 1))) ||u^r|| \
&= sqrt(((lambda^2 + 1 - 2 lambda^2) / (lambda^2 + 1)) / ((2 lambda^2 + 2 - 2 lambda^2) / (lambda^2 + 1))) ||u^r|| = sqrt((1 - (2 lambda^2) / (lambda^2 + 1)) / (2 - (2 lambda^2) / (lambda^2 + 1))) ||u^r|| = sqrt((1 - rho) / (2 - rho)) ||u^r||.
$ <eq:cos>

This idea can be generalized to higher-dimensional function spaces for the following reason: Since $u^(q - 1)$ is projection of $u^r$ onto $S_H$, the length of $u'$ is the shortest distance between $u^r$ and $S_H$. If we add more dimensions to $S_H$, the shortest distance can only become smaller but not larger. Therefore, @eq:cos still holds.

This equation gives us an upper bound for the norm reduction achieved by an exact coarse grid correction. Together with @eq:presmoother-convergence, we find
$
||u'|| &<= rho^(r / 2) sqrt((1 - rho) / (2 - rho)) ||u^0||
$ <eq:exact-cgc-convergence>
to account for both the presmoother and the exact coarse grid correction.

#{
let xs = lq.linspace(0, 1, num: 250)
let r = 3
show lq.selector(lq.legend): set grid(row-gutter: 1em)
place($(r = 3)$, dx: 107pt, dy: 11pt)
figure(
	lq.diagram(
		width: 200pt,
		xaxis: (
			label: $rho$
		),
		legend: (
			position: (100% + 0.6em, 0%)
		),
		lq.plot(
			xs,
			x => calc.pow(x, r / 2),
			label: $rho^(r/2)$,
			mark: none,
			stroke: 2pt,
		),
		lq.plot(
			xs,
			x => calc.sqrt((1 - x) / (2 - x)),
			mark: none,
			label: $sqrt((1 - rho) / (2 - rho))$,
			stroke: 2pt,
		),
		lq.plot(
			xs,
			x => calc.pow(x, r / 2) * calc.sqrt((1 - x) / (2 - x)),
			mark: none,
			label: $rho^(r/2) sqrt((1 - rho) / (2 - rho))$,
			stroke: 2pt,
		),
	),
	caption: [Norm reduction (red) after the presmoother (blue) and exact coarse grid correction (yellow) as a function of $rho$.]
)
}

Plotting the norm reduction of the presmoother and the exact coarse grid correction, we observe the following: In case the presmoother is not effective ($rho^(r / 2) = rho = 1$), the coarse grid correction is very effective ($sqrt((1 - rho) / (2 - rho)) -> 0$). Conversely, if the coarse grid correction is rather ineffective ($rho = 0$), the presmoother is very effective. The maximum of the red function, i.e., the worst possible norm reduction, is still far below $1$.

= Non-Exact Coarse Grid Correction <sec:non-exact-cgc>

Now we generalize the coarse grid correction to include an error. Instead of evaluating $u^(q - 1)$ exactly, we compute an approximation $v_1$ to $u^(q - 1)$. Assuming that we already know the convergence rate $delta_(q - 1)$ on the coarser level $q - 1$, and that perform a V-cycle with one recursive multigrid call, we find that the deviation from the exact error, $u^(q - 1)$ can be bounded:
$
||v_1 - u^(q - 1)|| <= delta_(q - 1) ||u^(q - 1)||.
$
To generalize this, we define
$
delta := cases(0 &wide& q = 1, (delta_(q - 1))^mu &wide& q > 1) wide wide wide mu := cases(1 &wide& "V-cycle", 2 &wide& "W-cycle")
$
and require that
$
||v_1 - u^(q - 1)|| <= delta ||u^(q - 1)||.
$ <eq:step-3-generalization>
Remember that $u^(r + 1)$ is the actual approximation after the coarse grid correction whereas $u'$ is the approximation if the coarse grid correction had been exact.
Thus the deviation is,
$
u^(r + 1) - u' = (u^r - v_1) - (u^r - u^(q - 1)) = u^(q - 1) - v_1.
$
If we define
$
v_2 := 1/delta (u^(q - 1) - v_1) = 1/delta (u^(r + 1) - u') wide wide v_2in S_H
$ <eq:v2>
then (from @eq:step-3-generalization):
$
||u^(r + 1) - u'|| = ||u^(q - 1) - v_1|| = delta ||v_2|| <= delta ||u^(q - 1)|| \
==> ||v_2|| <= ||u^(q - 1)||.
$ <eq:length-v2>

#figure(
	cetz.canvas({
		import cetz.draw: *

		scale(0.8)
		
		line((-3, 0), (8, 0))

		line((0, 0), (6, 6), stroke: 2pt, mark: (end: ">"))
		line((6, 0), (6, 6), stroke: 2pt, mark: (end: ">"))
		line((0, 0), (6, 0), stroke: 2pt, mark: (end: ">"))
		line((0, 0), (1.5, 0), stroke: 3pt, mark: (end: ">"))

		line((6, 6), (7.5, 6), stroke: (dash: "dashed"))
		line((6, 0), (7.5, 6), stroke: 2pt, mark: (end: ">"))
		line((1.5, 0), (7.5, 6), stroke: (dash: "dashed"))

		content((3, 3), $u^r$, anchor: "south-east")
		content((6, 0), $u^(q - 1)$, anchor: "north", padding: 0.3)
		content((0, 0), $0$, anchor: "north", padding: 0.2)
		content((6, 3), $u'$, anchor: "east", padding: 0.1)
		content((6.75, 3), $u^(r + 1)$, anchor: "west", padding: 0.2)
		content((8, 0), $S_H$, anchor: "west", padding: 0.2)
		content((1.5, 0), $delta v_2$, anchor: "north", padding: 0.2)

		right-angle((6, 0), start: 90deg)

		let h = -1.1
		let o = 0.15
		line((1.5, h), (6, h))
		line((1.5, h + o), (1.5, h - o))
		line((6, h + o), (6, h - o))
		content((3.75, h), $||v_1||$, anchor: "north", padding: 0.1)
	}),
	caption: [Visualization of the non-exact coarse grid correction.]
) <fig:non-exact-cgc>

= Full Multigrid Iteration <sec:full-mg>

To account for a full multigrid iteration, we need to include the presmoother, non-exact coarse grid correction and the postsmoother. First, let us analyze the algorithm across the first two iterations.

#figure(
	table(
		columns: 7,
		align: (right + horizon,) + (center + horizon,) * 6,
		stroke: none,
		[Iteration], table.vline(stroke: 2pt), table.cell(colspan: 3)[1], table.vline(stroke: 2pt), table.cell(colspan: 3)[2],
		table.hline(stroke: 2pt),
		[Step], [presmoother], [CGC], [postsmoother], [presmoother], [CGC], [postsmoother],
		table.hline(),
		[], $G^(r + 1) G^r ... G^1$, [], $G^1 ... G^r G^(r + 1)$, $G^r ... G^1$, [], $G^1...G^r G^(r + 1)$,
		table.hline(),
		[Example \ $r = 3$], $G2(G^4) G1(G^3) G2(G^2) G1(G^1)$, [], $G1(G^1) G2(G^2) G1(G^3) G2(G^4)$, $G1(G^3) G2(G^2) G1(G^1)$, [], $G1(G^1) G2(G^2) G1(G^3) G2(G^4)$,
	),
	caption: [First two iterations of the multigrid algorithm. \ The coarse grid correction is abbreviated as CGC.]
) <tab:no-duplicate>

Note that in the first iteration, we have an addition $G^(r + 1)$ at the beginning of the presmoother. This half-step is missing in all subsequent iterations. If we were to include $G^(r + 1)$ at the beginning of the presmoother of each iteration, we would have a duplicate. For example for $r = 3$ the postsmoother of iteration 1 would end with $G_h^"II"$ and the presmoother of iteration 2 would start with $G_h^"II"$ again. However, since the Gauss-Seidel half-steps are idempotent (as they are projections), we have $G_h^"II" dot G_h^"II" = G_h^"II"$. The additional half-step in the presmoother would not have any effect.

However, for the same reason, it is possible to include $G^(r + 1)$ conceptionally at the beginning of each presmoother. Mathematically, this is the same algorithm as in @tab:no-duplicate. This has the benefit that now the postsmoother is the reverse of the presmoother, such that
$
tilde(G) &= G^(r + 1) G^r ... G^1 &wide& "presmoother" \
tilde(G)^* &= G^1 ... G^r G^(r + 1) &wide& "postsmoother".
$ <eq:presmoother-postsmoother>
As suggested by the notation, $tilde(G)^*$ is the adjoint of $tilde(G)$:

$
(tilde(G) u, v) &= (G^(r + 1) G^r ... G^1 u, v) \
&= (G^r G^(r - 1) ... G^1 u, G^(r + 1) v) \
&= (G^(r - 1) ... G^1 u, G^r G^(r + 1) v) \
&= ... = (u, G^1 G^2 ... G^r G^(r + 1) v) \
&= (u, tilde(G)^* v).
$ <eq:G-adjoint>
This is true because all Gauss-Seidel half-steps $G^k$ are self-adjoint as shown in @sec:presmoother.

Now, we are ready to take the final steps to get the convergence rate. To this end, we use the following idea: We want to get an expression for $||u^(2 r + 2)||$. This term can be expressed as follows:
$
max_(hat(u) in S_h) {((hat(u), u^(2 r + 2))) / (||hat(u)||)} = max_(hat(u) in S_h) {(hat(u) / (||hat(u)||), u^(2 r + 2))} = ||u^(2 r + 2)||.
$ <eq:duality-start>
This is because (due to the Cauchy-Schwarz inequality, @eq:cauchy-schwarz)
$
(hat(u) / (||hat(u)||), u^(2 r + 2)) <= lr(||hat(u) / (||hat(u)||)||) dot ||u^(2 r + 2)|| = ||u^(2 r + 2)||.
$ <eq:duality>
Subsequently, we will find the upper bound
$
(hat(u) / (||hat(u)||), u^(2 r + 2)) <= delta_q ||u^0|| wide forall hat(u) in S_h.
$
Since this upper bound must be valid for all possible $hat(u)$, it must be valid vor $hat(u) = u^(2 r + 2)$ as well (when we have an equality in @eq:duality). In this special case, we find
$
(u^(2 r + 2) / (||u^(2 r + 2)||), u^(2 r + 2)) = ||u^(2 r + 2)|| <= delta_q ||u^0||.
$ <eq:duality-end>
To evaluate $delta_q$, we start with the inner product $(hat(u), u^(2 r + 2))$.

$
(hat(u), u^(2 r + 2)) &= (hat(u), tilde(G)^* u^(r + 1)) &wide& #ref(<eq:Gstar>) \
&= (tilde(G) hat(u), u^(r + 1)) &wide& #ref(<eq:G-adjoint>) \
&= (tilde(G) hat(u), u' + delta v_2) &wide& #ref(<eq:v2>) \
&= (tilde(G) hat(u), u' - delta u' + delta u' + delta v_2) \
&= (tilde(G) hat(u), (1 - delta) Q u^r + delta (u' + v_2)).
$
Here, $Q := P_(S_H)^perp$ represents the exact coarse grid correction (see @eq:exact-cgc). Since it is a projection, it is idempotent ($Q^2 = Q$).
$
(hat(u), u^(2 r + 2)) &<= (1 - delta) (tilde(G) hat(u), Q u^r) + delta (tilde(G) hat(u), u' + v_2) &wide& "linearity of inner product" \
&= (1 - delta) (tilde(G) hat(u), Q tilde(G) u^0) + delta (tilde(G) hat(u), u' + v_2) &wide& #ref(<eq:presmoother-postsmoother>) \
&<= (1 - delta) (tilde(G) hat(u), Q^2 tilde(G) u^0) + delta ||tilde(G) hat(u)|| ||u' + v_2|| &wide& #ref(<eq:cauchy-schwarz>).
$
For the term $||u' + v_2||$ we find (see @fig:non-exact-cgc and @eq:length-v2)
$
||u' + v_2||^2 = ||u'||^2 + ||v_2||^2 <= ||u'||^2 + ||u^(q - 1)||^2 = ||u^r||^2
$
and thus
$
||u' + v_2|| <= ||u^r||.
$
$
(hat(u), u^(2 r + 2)) &<= (1 - delta) (Q tilde(G) hat(u), Q tilde(G) u^0) + delta ||tilde(G) hat(u)|| ||u^r|| &wide& #ref(<eq:orthogonal-projections-self-adjoint>) \
&<= (1 - delta) ||Q tilde(G) hat(u)|| ||Q tilde(G) u^0|| + delta ||tilde(G) hat(u)|| ||tilde(G) u^0|| &wide& #ref(<eq:cauchy-schwarz>).
$

We rearrange this term as follows. If
$
a &:= sqrt(1 - delta) ||Q tilde(G) hat(u)||, &wide& c := sqrt(1 - delta) ||Q tilde(G) u^0||, \
b &:= sqrt(delta) ||tilde(G) hat(u)||, &wide& d := sqrt(delta) ||tilde(G) u^0||,
$
then
$
(1 - delta) ||Q tilde(G) hat(u)|| ||Q tilde(G) u^0|| + delta ||tilde(G) hat(u)|| ||tilde(G) u^0|| = a c + b d = vec(a, b) dot vec(c, d) \
<= lr(||vec(a, b)||) dot lr(||vec(c, d)||) = sqrt(a^2 + b^2) dot sqrt(c^2 + d^2) \
= sqrt((1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2) dot sqrt((1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2).
$ <eq:sqrt-terms>
Let us focus on the term $(1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2$. From @eq:exact-cgc-convergence and @eq:presmoother-convergence, we find
$
(1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2 &= (1 - delta) ||u'||^2 + delta ||u^r||^2 \
&<= (1 - delta) (rho^(r/2) sqrt((1 - rho) / (2 - rho)) ||u^0||)^2 + delta (rho^(r/2) ||u^0||)^2 \
&= (1 - delta) rho^r (1 - rho) / (2 - rho) ||u^0||^2 + delta rho^r ||u^0||^2 \
&= rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta) ||u^0||^2.
$ <eq:rho-u0>
Similarly, for the other square root in @eq:sqrt-terms, we have
$
(1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2 <= rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta) ||hat(u)||^2.
$ <eq:rho-hatu>

However, $rho$ depends on the starting vector (previously called $u^0$). When we have different starting vectors, $u^0$ and $hat(u)$, the corresponding $rho$ will be different (i.e. the value of $rho$ in @eq:rho-u0 may be different than in @eq:rho-hatu). Furthermore, we cannot evaluate the exact value of $rho$ since we can start with an arbitrary start vector $u^0$. To overcome this problem, we assume the worst case, i.e., the maximum of the function $rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)$.

To do so, we find that $0 <= rho <= 1$. This is because
$
0 <= lambda := (|v^r|) / (||v^r||) <= 1 wide "and" wide rho := (2 lambda^2) / (lambda^2 + 1).
$


#figure({
	let xs = lq.linspace(0, 3, num: 200)
	lq.diagram(
		xaxis: (
			label: $lambda$,
		),
		yaxis: (
			label: $rho$,
		),
		lq.plot(
			xs,
			x => 2 * x * x / (x * x + 1),
			mark: none,
		)
	)
}, caption: [$rho$ as a function of $lambda$.])

We define
$
delta_q := max_(0 <= rho <= 1) [rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)].
$
Then
$
(1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2 &<= delta_q ||u^0||^2 \
(1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2 &<= delta_q ||hat(u)||^2.
$
We substitute this inequality back into @eq:sqrt-terms:
$
(hat(u), u^(2 r + 2)) &<= sqrt(delta_q ||hat(u)||^2) dot sqrt(delta_q ||u^0||^2) \
&= delta_q ||hat(u)|| ||u^0|| wide wide forall hat(u) in S_h.
$
Using the arguments from @eq:duality-start to @eq:duality-end, we find that
$
||u^(2 r + 2)|| <= delta_q ||u^0|| wide "with" wide delta_q := max_(0 <= rho <= 1) [rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)].
$
$delta_q$ is an upper bound for the norm reduction of a single multigrid iteration on level $q$. This equation can now be applied for both a two-grid method and a general multigrid method to get concrete numerical values for $delta_q$.
