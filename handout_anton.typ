#import "lib.typ": *

#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#set par(justify: true)
#set page(numbering: "1")

#let hi(content) = box(width: 100%, inset: 0.8em, stroke: orange, content)
#let infobox(content) = box(width: 100%, inset: 0.8em, stroke: blue, content)

#title[Convergence Rate of a Multigrid Method with Gauss-Seidel Relaxation for the Poisson Equation]

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

#let draw-2h-grid() = {
	import cetz.draw: line
	for x in range(0, W + 1, step: 2) {
		line((x, 0), (x, H))
	}
	for y in range(0, H + 1, step: 2) {
		line((0, y), (W, y))
	}
	for x in range(0, W, step: 2) {
		for y in range(0, H, step: 2) {
			if odd(int(x / 2) + int(y / 2)) {
				line((x + 2, y), (x, y + 2))
			} else {
				line((x, y), (x + 2, y + 2))
			}
		}
	}
	for x in range(0, W + 1, step: 2) {
		for y in range(0, H + 1, step: 2) {
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


= Setup

Poisson equation with Dirichlet boundary conditions
$
-laplace u = f &wide& &"in" Omega \
u = 0 &wide& &"on" partial Omega
$

= Grids

$
Omega_h &: "Dreiecke mit Kantenlängen" h "und" sqrt(2) h \
Omega_H &: "Dreiecke mit Kantenlängen" sqrt(2) h "und" 2 h \
Omega_(2h) &: "Dreiecke mit Kantenlängen" 2h "und" 2 sqrt(2) h
$

#figure(
	cetz.canvas({
		import cetz.draw: content, translate, scale
		scale(0.7)
		draw-grid()
		content((W / 2, H + 0.5), $Omega_h$)
		translate(x: 8)
		draw-coarse-grid()
		content((W / 2, H + 0.5), $Omega_H$)
		translate(x: 8)
		draw-2h-grid()
		content((W / 2, H + 0.5), $Omega_(2h)$)
	}),
	caption:["The grid is divided into fine points (blue) and coarse points (yellow)"]
)

= Decomposition

$
S_h &: "finite Elemente Raum auf" Omega_h \
S_H &: "finite Elemente Raum auf" Omega_H \
T_h &:= {w in S_h | w(p) = 0 wide forall p in Omega_H}
$

#figure(
	cetz.canvas({
		import cetz.draw: content, translate, scale
		scale(0.7)
		draw-grid()
		content((7, 2), $=$)
		content((W/2, -0.5), $u in S_h$)
		translate(x: 8)
		draw-coarse-grid()
		content((7, 2), $plus.o$)
		content((W/2, -0.5), $v in S_H$)
		translate(x: 8)
		draw-th-grid()
		content((W/2, -0.5), $w in T_h$)
	}),
	caption: [The space $S_h$ of the fine grid can be expressed as a sum of the space $S_H$ on the coarse grid (middle) and the space $T_h$ (right).]
)

$
S_h &= S_H plus.o T_h wide u &= v plus w
$

The basis functions of $T_h$ form pyramid-like patterns. As a result, all functions in $T_h$ are zero on the diagonal lines through the coarse grid points.

#figure(
	cetz.canvas({
		cetz.draw.scale(0.8)
		draw-grid(show-th: true)
	}),
	caption: [Each shaded are corresponds to one of the basis functions of $T_h$. The supports of the respective basis functions are disjoint (they do not overlap each other).]
)
= Norms

For $u, v in S_h$, we define the bilinear form
$
a(u, v) := integral_Omega (u_xi v_xi + u_eta v_eta) thin dif xi thin dif eta.
$

Furthermore, the energy norm is given by
$
||v|| := sqrt(a(v, v)).
$

Therefore,
$
||v||^2 = a(v, v) &= integral_Omega v_xi^2 + v_eta^2 \
&= sum_nu integral_nu v_xi^2 + v_eta^2 \
&= sum_nu (lr(v_xi^2|)_nu + lr(v_eta^2|)_nu) underbrace(1/2 h^2, op("area")(nu)) \
&= sum_nu ((v_"east" - v_"west")^2 / h^2 + (v_"north" - v_"south")^2 / h^2) 1/2 h^2 \
&= sum_nu 1/2 (v_"east" - v_"west")^2 + sum_nu 1/2 (v_"north" - v_"south")^2 \
&= sum_vec(i\, j in Omega_h, d(i, j) = h, delim: #none) (v_i - v_j)^2.
$


The last step is justified because each vertical and horizontal border is adjacent to two different triangles (accounting for the factor 2). Note that each pair $(i, j)$ is only used once, meaning that if we include pair $(i, j)$ in the sum ($i != j$), we do not include $(j, i)$ as well.

Similarly, we define a seminorm (i.e., something that looks like a norm but that could evaluate to zero, $|x| = 0$, even though $x != 0$),

$
|v|^2 := sum_vec(i\, j in Omega_H, d(i, j) = 2 h, delim: #none) 1/2 (v_i - v_j)^2.
$

#hi[Note that for both norms the sum is only taken over all vertical and horizontal edges of the triangulation (i.e., pairs $(i, j)$ with Euclidean distance $h$), not over the diagonal edges.]

For any $p_i, p_j in Omega_H$ with $d(i, j) = 2 h$, the midpoint $p_k$ between them lies in $Omega_h$ and satisfies $d(i, k) = d(k, j) = h$. Since
$
(u_i - u_j)^2 = ((u_i - u_k) + (u_k - u_j))^2 <= 2(u_i - u_k)^2 + 2(u_k - u_j)^2
$
(this follows from $((u_i - u_k) - (u_k - u_j))^2 >= 0$), summing over all such pairs gives
$
|u|^2 = sum_vec(i\, j in Omega_H, d(i, j) = 2 h, delim: #none) 1/2 (u_i - u_j)^2 <= sum_vec(i\, j in Omega_H, d(i, j) = 2 h, delim: #none) (u_i - u_k)^2 + (u_k - u_j)^2.
$
Every point $p_k in Omega_h \\ Omega_H$ has exactly two $2 h$-edges of $Omega_H$ passing through it, contributing exactly its four incident $h$-edges of $Omega_h$. Hence, summing over all pairs $(i, j) in Omega_H$, every $h$-edge $(i, k)$ with $p_k in Omega_h \\ Omega_H$ appears exactly once on the right-hand side, so it equals $||u||^2$. Therefore,
$
|u| <= ||u|| wide forall u in S_h.
$

= Estimating a(v,w)

#figure(
	cetz.canvas({
		import cetz.draw: line, content, scale, circle

		scale(2)

		let p1 = (0, -1)
		let p2 = (1, 0)
		let p3 = (0, 1)
		let p4 = (-1, 0)
		let p5 = (0, 0)

		line((p4.at(0) - 0.5, 0), (p2.at(0) + 0.5, 0), stroke: (dash: "dashed", thickness: 0.5pt), mark: (end: ">"), name: "xi-axis")
		line((0, p1.at(1) - 0.3), (0, p3.at(1) + 0.4), stroke: (dash: "dashed", thickness: 0.5pt), mark: (end: ">"), name: "eta-axis")
		content("xi-axis.end", $xi$, anchor: "north-west", padding: 0.1)
		content("eta-axis.end", $eta$, anchor: "south-east", padding: 0.1)

		line(p1, p2, p3, p4, close: true)
		line(p4, p2)

		let r = 2pt
		circle(p1, radius: r, stroke: none, fill: red)
		circle(p2, radius: r, stroke: none, fill: olive)
		circle(p3, radius: r, stroke: none, fill: red)
		circle(p4, radius: r, stroke: none, fill: olive)
		circle(p5, radius: r, stroke: none, fill: black)

		content(p1, $1$, anchor: "east", padding: 0.25)
		content(p3, $3$, anchor: "east", padding: 0.25)
		content(p2, $2$, anchor: "north", padding: 0.25)
		content(p4, $4$, anchor: "north", padding: 0.25)
		content(p5, $5$, anchor: "north-east", padding: 0.2)

		content((0.5, 0), $h$, anchor: "south", padding: 0.1)
		content((0, 0.5), $h$, anchor: "west", padding: 0.1)
	})
)

Here, $"I" = triangle 1 2 4$ and $"II" = triangle 4 2 3$.

$
a(v,w)_"I" = integral_"I" v_xi w_xi + v_eta w_eta = integral_"I" v_xi w_xi + integral_"I" v_eta w_eta
$

$
integral_"I" v_xi w_xi = integral_(I_"links") v_xi w_xi + integral_(I_"rechts") v_xi w_xi = 0
$
da $lr(w_xi|)_(I_"links") = - lr(w_xi|)_(I_"rechts")$ und $v_xi = "const"$.

$
integral_"I" v_eta w_eta
&= integral_"I" 1/2 v_eta^2 - 1/2 v_eta^2 + 1/2 w_eta^2 - 1/2 w_eta^2 + v_eta w_eta \
&= integral_"I" 1/2 (-v_eta^2 + 2 v_eta w_eta - w_eta^2) + 1/2 v_eta^2 + 1/2 w_eta^2 \
&= integral_"I" 1/2 v_eta^2 + 1/2 w_eta^2 - 1/2 (v_eta^2 - 2 v_eta w_eta + w_eta^2) \
&= 1/2 integral_"I" v_eta^2 + 1/2 integral_"I" w_eta^2 - 1/2 integral_"I" (v_eta + w_eta)^2 \
&= 1/2 integral_"I" v_eta^2 + 1/2 integral_"I" v_xi^2 - 1/2 integral_"I" v_xi^2 + 1/4 integral_"I" (w_eta^2 + w_eta^2) - 1/2 integral_"I" (v_eta + w_eta)^2 \
&= underbrace(1/2 integral_"I" (v_eta^2 + v_xi^2), A_"I") + underbrace(1/4 integral_"I" (w_eta^2 + w_xi^2), B_"I") - underbrace(1/2 integral_"I" v_xi^2, C_"I") - underbrace(1/2 integral_"I" (v_eta + w_eta)^2, D_"I") \
&= A_"I" + B_"I" - C_"I" - D_"I"
$

When summing up the term $A$ for all triangles in the triangulation, we simply get
$
1/2 integral_Omega v_eta^2 + v_xi^2 =1/2 thin a(v, v) = 1 / 2 thin ||v||^2.
$ <eq:A>

Similarly, for $B$:
$
1/4 integral_Omega w_eta^2 + w_xi^2 = 1/4 thin a(w, w) = 1/4 thin ||w||^2.
$ <eq:B>

Since $v$ is linear in both $"I" = triangle 1 2 4$ and in $"II" = triangle 4 2 3$, $v_xi$ is a constant in both triangles. However, because of the shared boundary $overline(4 2)$, we know that $lr(v_xi|)_"I" = lr(v_xi|)_"II" = (v_2 - v_4) / (2 h)$.

For $C$ we find
$
sum_nu C_nu = sum_vec(i\, j "green", d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2
$

We can also show that
$
sum_nu D_nu >= sum_vec(i\, j "red", d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2
$

Together with $C$, we get
$
sum_nu C_nu + D_nu &>= sum_vec(i\, j "green", d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2 + sum_vec(i\, j "red", d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2 \
&= sum_vec(i\, j in Omega_H, d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2 \
&= 1/2 sum_vec(i\, j in Omega_H, d(i, j) = 2 h, delim: #none) 1/2 (v_i - v_j)^2 \
&= 1/2 |v|^2
$ <eq:CD>
according to the definition of the seminorm $|dot|$.


When summing up the terms for $A$ (@eq:A), $B$ (@eq:B), $C$ and $D$ (@eq:CD), we find that

$
a(v, w) <= underbrace(1/2 ||v||^2, A) + underbrace(1/4 ||w||^2, B) - underbrace(1/2 |v|^2, C + D) = 1/2 (||v||^2 - |v|^2) + 1/4 ||w||^2.
$

#hi[Notice that we have an inequality ($<=$) here because of the bound on $sum_nu D_nu$.]

Interestingly, this upper bound is quadratic in $||w||$ but at the same time, bilinear form $a(v, w)$ is linear in $w$ since $a(v, t w) = t thin a(v, w)$ for $t in RR$ (homogeneity). If we set $w =: t hat(w)$ such that $||hat(w)|| = 1$ and $t > 0$, we get
$
a(v, t hat(w)) =  t thin a(v, hat(w)) &<= 1/2 (||v||^2 - |v|^2) + 1/4 t^2 wide | : t \
a(v, hat(w)) &<= (||v||^2 - |v|^2) / (2 t) + 1/4 t.
$

We substitute $mu := ||v||^2 - |v|^2$ for brevity. This function has the following minimum
$
pdv(, t) (mu / (2 t) + 1/4 t) = -mu / (2 t^2) + 1/4 =^! 0 \
t_* = sqrt(2 mu).
$

This is the value where the tangent $f(t) = C t$ (yellow) touches the quadratic function $1/2 mu + 1/4 t^2$ (blue). The derivative is shown in red.

#figure(
	lq.diagram(
		yaxis: (
			position: 0,
			lim: (-3, 5),
			ticks: none,
			tip: tiptoe.triangle,
		),
		xaxis: (
			position: 0,
			ticks: none,
			extra-ticks: (lq.tick(1, label: $t_*$, inset: 2pt, outset: 2pt),),
			tip: tiptoe.triangle,
			label: lq.label($t$, dx: 85pt, dy: -16pt),

		),
		lq.plot(
			lq.linspace(-5, 5),
			x => 1 + x * x,
			mark: none,
		),
		lq.plot(
			lq.linspace(-5, 0, include-end: false),
			x => 1/x + x,
			mark: none,
			stroke: lq.color.map.petroff10.at(2)
		),
		lq.plot(
			lq.linspace(0.1, 5),
			x => 1/x + x,
			mark: none,
			stroke: lq.color.map.petroff10.at(2)
		),
		lq.plot(
			lq.linspace(-5, 5),
			x => 2 * x,
			mark: none,
			stroke: lq.color.map.petroff10.at(1)
		),
		lq.line(
			(1, 0),
			(1, 2),
			stroke: (dash: "dashed", paint: gray),
		),
	)
)


Accordingly,
$
a(v, hat(w)) &<= mu / (2 t_*) + 1/4 t_* \
&= mu / (2 sqrt(2 mu)) + 1/4 sqrt(2 mu) \
&= 1/sqrt(8) sqrt(mu) + sqrt(2) / sqrt(16) sqrt(mu) \
&= sqrt(1/2 mu).
$

Now we multiply on both sides by $t = ||w|| > 0$ and we get
$
t thin a(v, hat(w)) = a(v, t thin hat(w)) = a(v, w) <= sqrt(1/2 mu) ||w||.
$

Because of homogeneity, when we negate $w$, we get
$
-a(v, w) = a(v, -w) <= sqrt(1/2 mu) ||-w|| = sqrt(1/2 mu) ||w||.
$
Therefore, we can replace $a(v, w)$ by its absolute value. Finally, this yields
$
|a(v, w)| <= sqrt(1/2 thin (||v||^2 - |v|^2)) ||w|| = sqrt(1/2 thin (1 - (|v|^2) / (||v||^2))) ||v||||w||.
$

*Geometric Interpretation*


#figure(
	cetz.canvas({
		import cetz.draw: line, content, arc
		let a = (0, 0)
		let b = (6, 0)
		let c = (4, 2)
		line(a, b, mark: (end: ">"), name: "v")
		line(b, c, mark: (end: ">"), name: "w")
		line(a, c, mark: (end: ">"), name: "u")
		content("v", $v$, padding: 0.1, anchor: "north")
		content("u", $u$, padding: 0.1, anchor: "south-east")
		content("w", $w$, padding: 0.1, anchor: "south-west")
		arc(b, anchor: "origin", radius: 1, start: 135deg, stop: 180deg)
		content((5.4, 0.28), $alpha$)
	})
)

We decompose $u$ into the sum $v + w$. These vectors are part of the same Hilbert space with inner product $a(dot, dot)$. We define the angle $alpha$ as
$
cos(alpha) = a(v, w)/(||v||||w||).
$
Therefore, lemma 3.1 can be written in the following way,
$
|cos(alpha)| = (|a(v, w)|) / (||v||||w||) <= sqrt(1/2 (1 - ((|v|) / (||v||))^2))
$
or, alternatively,
$
cos^2 (alpha) <= 1/2 (1 - lambda^2)
$
where $lambda := (|v|) / (||v||)$. Due to the identity $sin^2 alpha + cos^2 alpha = 1$, we can rearrange this inequality,
$
cos^2 alpha = 1 - sin^2 alpha <= 1/2 (1 - lambda^2) wide ==> wide sin^2 alpha >= 1/2 (1 + lambda^2).
$
Furthermore, we observe that for a fixed $v$ and a fixed $alpha$, $u$ is the shortest when $u perp w$. In this case, $sin(alpha) = (||u||) / (||v||)$. Generally, $||u|| >= ||v|| sin(alpha)$.
Since $u$ and $v$ have the same values on the green and red points, we know that $|u| = |v|$.
$
|u| = |v| = lambda||v|| <= lambda / sin(alpha)||u|| <= lambda / (sqrt(1/2 (1 + lambda^2)))||u||
$
Finally, we express this as
$
(|u|) / (||u||) <= lambda / sqrt(1/2 (1 + lambda^2)) = sqrt(lambda^2 / (1/2 (1 + lambda^2))) = sqrt((2 lambda^2) / (1 + lambda^2)) =: sqrt(rho)
$

= Gauß-Seidel Relaxation

We use a finite element discretization. Interestingly, the stencil for the black all points is the same, wheter in $Omega_h$ or in $Omega_H$:

$
mat(0, -1, 0; -1, 4, -1; 0, -1, 0)
$



#figure(
cetz.canvas({
	import cetz.draw: line, circle, translate, content, floating

	let p0 = ( 0,  0)
	let p1 = (-1, -1)
	let p2 = ( 0, -1)
	let p3 = ( 1, -1)
	let p4 = ( 1,  0)
	let p5 = ( 1,  1)
	let p6 = ( 0,  1)
	let p7 = (-1,  1)
	let p8 = (-1,  0)

	let d(coords, color) = circle(coords, fill: color, stroke: none, radius: 2pt)
	let coordinate-system() = {
		line((-1.7, 0), (1.7, 0), stroke: (dash: "dashed", thickness: 0.5pt), mark: (end: ">"), name: "xi")
		floating(content("xi.end", $xi$, anchor: "north-west", padding: 0.05))
		line((0, -1.7), (0, 1.7), stroke: (dash: "dashed", thickness: 0.5pt), mark: (end: ">"), name: "eta")
		floating(content("eta.end", $eta$, anchor: "east", padding: 0.11))
	}
	let point-labels() = {
		content(p2, $1$, anchor: "north-west", padding: 0.05)
		content(p4, $2$, anchor: "north-west", padding: 0.05)
		content(p6, $3$, anchor: "south-east", padding: 0.05)
		content(p8, $4$, anchor: "south-east", padding: 0.05)
	}
	

	coordinate-system()
	line(p1, p3, p5, p7, close: true)
	line(p4, p8)
	line(p2, p6)
	line(p2, p4, p6, p8, close: true)
	d(p1, black)
	d(p2, olive)
	d(p3, black)
	d(p4, red)
	d(p5, black)
	d(p6, olive)
	d(p7, black)
	d(p8, red)
	d(p0, black)
	point-labels()
	content(p0, $0$, anchor: "north-west", padding: 0.05)
	
	translate(x: 6)
	
	coordinate-system()
	line(p1, p3, p5, p7, close: true)
	line(p3, p7)
	line(p1, p5)
	line(p4, p8)
	line(p2, p6)
	d(p1, olive)
	d(p2, black)
	d(p3, olive)
	d(p4, black)
	d(p5, olive)
	d(p6, black)
	d(p7, olive)
	d(p8, black)
	d(p0, red)
	point-labels()
	content((0.27, 0.37), $0$, anchor: "north-west", padding: 0.05)
}),
caption: [Neighborhood of a fine grid node in $Omega_h without Omega_H$ (left) and of a coarse node $Omega_H$ (right).]
) <fig:neighborhood>

The stencil results in the following Gauss-Seidel relaxation for a single row:

$
u_i = 1/4 sum_(j = 1)^4 u_j + b_i
$
where $b_i = integral_Omega f phi_i$ is the right-hand side. In the following section, we will be concerned with the error of our approximation for $u$. Therefore, we assume that $b = 0$.

We split the Gauss-Seidel relaxation into two steps ($eq.est$ red-black Gauss-Seidel). First, we will smooth the fine points ($in Omega_h without Omega_H$), then the coarse points ($in Omega_H$). Formally,
$
(G_h^"I" u)_i &= cases(u_i &wide& p_i in Omega_H &wide& , 1/4 sum_(j in N(i)) u_j &wide& p_i in Omega_h \\ Omega_H &wide& 
) \
(G_h^"II" u)_i &= cases(1/4 sum_(j in N(i)) u_j &wide& p_i in Omega_H &wide& , u_i &wide& p_i in Omega_h \\ Omega_H &wide& )
$

Here, $N(i)$ is the set of neighboring points with distance h to $p_i$.

*Lemma 4.1*

Assume that $Omega$ is convex. We want to show: if $u = G_h^"I" u$ then
$
 ||G_h^"II" u|| <= |u|.
$

half a Gauss-Seidel step decreases the energy norm beyond the seminorm - this means it is reduced by at least the factor $sqrt(rho)$.

*Proof.* For convenience we first ignore the boundary $partial Omega$; we comment on the necessary modifications at the end of the proof.

Let $macron(u) := G_h^"II" u$, and fix a point $p_0 in Omega_H$ with its four neighbors $1,2,3,4$ at distance $h$, which necessarily lie in $Omega_h without Omega_H$ (cf. @fig:neighborhood, right).

#figure(
	cetz.canvas({
		import cetz.draw: line, circle, content, floating

		let p0 = ( 0,  0)
		let p1 = (-1, -1)
		let p2 = ( 0, -1)
		let p3 = ( 1, -1)
		let p4 = ( 1,  0)
		let p5 = ( 1,  1)
		let p6 = ( 0,  1)
		let p7 = (-1,  1)
		let p8 = (-1,  0)

		let d(coords, color) = circle(coords, fill: color, stroke: none, radius: 2pt)

		line((-1.7, 0), (1.7, 0), stroke: (dash: "dashed", thickness: 0.5pt), mark: (end: ">"), name: "xi")
		floating(content("xi.end", $xi$, anchor: "north-west", padding: 0.05))
		line((0, -1.7), (0, 1.7), stroke: (dash: "dashed", thickness: 0.5pt), mark: (end: ">"), name: "eta")
		floating(content("eta.end", $eta$, anchor: "east", padding: 0.11))

		line(p1, p3, p5, p7, close: true)
		line(p3, p7)
		line(p1, p5)
		line(p4, p8)
		line(p2, p6)
		d(p1, olive)
		d(p2, black)
		d(p3, olive)
		d(p4, black)
		d(p5, olive)
		d(p6, black)
		d(p7, olive)
		d(p8, black)
		d(p0, red)
		content(p2, $1$, anchor: "north-west", padding: 0.05)
		content(p4, $2$, anchor: "north-west", padding: 0.05)
		content(p6, $3$, anchor: "south-east", padding: 0.05)
		content(p8, $4$, anchor: "south-east", padding: 0.05)
		content((0.27, 0.37), $0$, anchor: "north-west", padding: 0.05)
	}),
	caption: [Point $p_0 in Omega_H$ and its four neighbors $1,2,3,4 in Omega_h without Omega_H$ at distance $h$]
) <fig:lemma41-diamond>

For arbitrary real numbers $z_1, z_2, z_3, z_4$ one has the algebraic identity
$
2 sum_(j=1)^4 z_j^2 = \ (z_1 - z_2)^2 + (z_2 - z_3)^2 + (z_3 - z_4)^2 + (z_4 - z_1)^2 + 1/2 (z_1 + z_2 + z_3 + z_4)^2 - 1/2 (z_1 - z_2 + z_3 - z_4)^2.
$ <eq:diamond-identity>

Since $macron(u) = G_h^"II" u$, the value at $p_0$ is by definition the mean of the (unaffected) values at its neighbors,
$
macron(u)_0 = 1/4 sum_(i=1)^4 macron(u)_i = 1/4 sum_(i=1)^4 u_i,
$
because $G_h^"II"$ leaves the values at $Omega_h without Omega_H$ unchanged. Put $z_i := macron(u)_i - macron(u)_0 = u_i - macron(u)_0$. Then
$
sum_(j=1)^4 z_j = sum_(j=1)^4 u_j - 4 macron(u)_0 = 0,
$
so the first-order term $1/2 (z_1+z_2+z_3+z_4)^2$ in @eq:diamond-identity vanishes, and dropping the remaining non-positive term $-1/2(z_1-z_2+z_3-z_4)^2$ only weakens the identity into an inequality. We obtain
$
sum_(j=1)^4 (macron(u)_j - macron(u)_0)^2 <= 1/2 [(u_1 - u_2)^2 + (u_2 - u_3)^2 + (u_3 - u_4)^2 + (u_4 - u_1)^2].
$ <eq:step1>

Recall from the definition of the energy norm (2.3) that, splitting the sum according to which grid the points belong to,
$
||macron(u)||^2 = sum_(i in Omega_H) sum_(j in Omega_h without Omega_H \ d(i,j) = h) (macron(u)_i - macron(u)_j)^2.
$
For fixed $i = p_0$, the inner sum is exactly the left-hand side of @eq:step1. Summing the bound @eq:step1 over all $p_0 in Omega_H$ therefore gives
$
||G_h^"II" u||^2 <= sum_(k,l in Omega_h without Omega_H \ d(k,l) = H) (u_k - u_l)^2,
$ <eq:step2>
where $H = sqrt(2) h$ is the distance between two opposite corners of the diamond in @fig:lemma41-diamond.

#hi[Every pair $(k,l) in Omega_h without Omega_H$ with $d(k,l) = H$ occurs in the diamonds of *two* neighboring points $p_0 in Omega_H$. The factor $1/2$ on the right-hand side of @eq:step1 exactly compensates for this double counting, so no extra factor appears in @eq:step2.]

This is the only place where the hypothesis $u = G_h^"I" u$ enters. Since $u = G_h^"I" u$, every point $k in Omega_h without Omega_H$ already satisfies the relaxation equation
$
u_k = 1/4 sum_(j in N(k)) u_j,
$
i.e. $u_k$ is the mean of its four neighbors, which lie in $Omega_H$.

Consider a pair $k, l in Omega_h without Omega_H$ with $d(k,l) = H$, and let $1,3,5,7 in Omega_H$ be enumerated so that $1, 5$ are the neighbors shared in the relevant direction (cf. Figure 4 of Braess' paper). Using the mean-value property at $k$ and $l$,
$
16 (u_k - u_l)^2 = (u_1 + u_5 - u_3 - u_7)^2 <= 2(u_1 - u_3)^2 + 2(u_5 - u_7)^2,
$
where the last step is the elementary inequality $(a+b)^2 <= 2a^2 + 2b^2$ applied to $a = u_1 - u_3$, $b = u_5 - u_7$.

Summing this over all pairs $(k,l)$ from @eq:step2 — each pair $(u_1 - u_3)^2$ now reappears for two different pairs $(k,l)$, which compensates the factor $1/4$ — we obtain
$
||G_h^"II" u||^2 <= 1/4 sum_(m,n in Omega_H \ d(m,n) = 2H) (u_m - u_n)^2.
$ <eq:step3>
The points $m,n$ are again in $Omega_H$, but now at distance $2H = 2 sqrt(2) h$: these are exactly the two diagonals of a square of side length $2h$ with corners in $Omega_H$.

It remains to bound such diagonals by sides of length $2h$. Let $1,2,3,4 in Omega_H$ be the corners of such a square, in cyclic order, so that $(u_1-u_3)^2$ and $(u_2-u_4)^2$ are the two diagonal terms appearing in @eq:step3. A direct expansion gives the identity
$
(u_1 - u_3)^2 + (u_2 - u_4)^2 = \ (u_1-u_2)^2 + (u_2-u_3)^2 + (u_3-u_4)^2 + (u_4-u_1)^2 - (u_1 - u_2 + u_3 - u_4)^2,
$
and dropping the (non-positive) last term yields
$
(u_1-u_3)^2 + (u_2-u_4)^2 <= (u_1-u_2)^2 + (u_2-u_3)^2 + (u_3-u_4)^2 + (u_4-u_1)^2,
$
i.e. the sum of the two diagonals is bounded by the sum of the four sides of the square (each of length $2h$).

Every side of length $2h$ is shared by exactly two neighboring squares, so summing this inequality over all squares with corners in $Omega_H$ and combining with @eq:step3 gives
$
||G_h^"II" u||^2 <= 1/2 sum_(n,m in Omega_H \ d(n,m) = 2h) (u_n - u_m)^2 = |u|^2,
$
which is exactly the definition of the seminorm $|dot|$ from (2.4). This proves
$
||G_h^"II" u|| <= |u|.
$

Near $partial Omega$, some of the neighbors $1,2,3,4$ (or $1,3,5,7$) used above may not exist inside $Omega$; by the homogeneous Dirichlet condition, the corresponding values are simply replaced by $0$. The same telescoping arguments go through verbatim with these zero values inserted — only the index bookkeeping becomes more tedious, while every inequality used above (@eq:diamond-identity, the mean-value property, the diagonal identity) remains valid for boundary configurations as well, consistent with how the doubled weight in (2.5) was introduced for the seminorm. $qed$

As a consequence of Lemma 4.1 we have the following corollary.

*Corollary 4.2.* If $Omega$ is convex, then $||G_h^"II" dot G_h^"I" u|| <= |u|$ for every $u in S_h$ — i.e. the hypothesis $u = G_h^"I" u$ in Lemma 4.1 can be dropped, at the price of applying $G_h^"I"$ first.

*Proof.* Put $u' := G_h^"I" u$. Since $G_h^"I"$ is idempotent ($G_h^"I" G_h^"I" = G_h^"I"$), we have $u' = G_h^"I" u'$, so Lemma 4.1 applies to $u'$:
$
||G_h^"II" u'|| <= |u'|.
$
Now $G_h^"I"$ only changes the $w$-part of the decomposition $u = v+w$ (it leaves the values at $Omega_H$, i.e. the $v$-part, untouched), and the seminorm $|dot|$ depends only on the $v$-part by (3.3). Hence
$
|u'| = |G_h^"I" u| = |u|,
$
and combining the two displays gives $||G_h^"II" G_h^"I" u|| <= |u|$ for every $u in S_h$. $qed$

#hi[We note that the constant $1$ in (4.2) is best possible, since the convergence rate of the (pure) Gauss-Seidel iteration is known to be close to $1$, with $||G_h^"II" dot G_h^"I"|| = 1 - O(h^2)$.]

= Alternate Method

Proof that $G^"I" = I - P_V$ where $P_V$ is an orthogonal projection and $V = T_h$. The first Gauß-Seidel half-step is associated with the splitting $S_h = S_H plus.o T_h$ where $T_h = { w in S_h | w(p_i) = 0 "for" p_i in Omega_H}$. All functions in $u + T_h$ have the same values as $u$ at the white nodes (and possibly different values at the black nodes).

#hi[Don't call it $w$ since $w in T_h$ and $G^"I" u in T_h^perp$.]

$w := G^"I" u$ minimizes $J(v)$ over $u + T_h$. This is because of the following reason:

$
J(v) = a(v, v) = ||v||^2 = sum_vec(i\,j, d(i, j) = h, delim: #none) (v_i - v_j)^2.
$
Now, if we look at the support of a single basis function of $T_h$ (see @fig:neighborhood left), we find that the relevant terms of the sum are
$
(v_1 - v_0)^2 + (v_2 - v_0)^2 + (v_3 - v_0)^2 + (v_4 - v_0)^2 = "function"(v_0),
$
and the minimum of this function is at $v_0 = 1/4 (v_1 + v_2 + v_3 + v_4)$, which is exactly what $G^"I"$ computes for $v_0$.

Since $w$ minimizes $J(v)$ over $u + T_h$, going in any direction $t in T_h$ must increase the energy: $J(w + t) >= J(w) space forall t in T_h$. This can also be written as
$
dv(, epsilon) J(w + epsilon t)|_(epsilon = 0) =^! 0 wide forall t in T_h.
$
Now we use the definition of $J(v)$ and $a(u, v)$ (we assume that the right-hand side $f = 0$).

#hi[no need to substitute the definition of $a$ since we know it is linear and symmetric]

$
dv(, epsilon) J(w + epsilon t)|_(epsilon = 0) &= dv(, epsilon) a(w + epsilon t, w + epsilon t)|_(epsilon = 0) \
&= dv(, epsilon) integral_Omega (pdv((w + epsilon t), xi))^2 + (pdv((w + epsilon t), eta))^2|_(epsilon = 0) \
&= dv(, epsilon) integral_Omega (w_xi + epsilon t_xi)^2 + (w_eta + epsilon t_eta)^2|_(epsilon = 0) \
&= dv(, epsilon) integral_Omega (w_xi^2 + 2 epsilon w_xi t_xi + epsilon^2 t_xi^2) + (w_eta^2 + 2 epsilon w_eta t_eta + epsilon^2 t_eta^2)|_(epsilon = 0) \
&= dv(, epsilon) [integral_Omega (w_xi^2 + w_eta^2) + 2 epsilon integral_Omega (w_xi t_xi + w_eta t_eta) + epsilon^2 integral_Omega (t_xi^2 + t_eta^2)]_(epsilon = 0) \
&= dv(, epsilon) [a(w, w) + 2 epsilon a(w, t) + epsilon^2 a(t, t)]_(epsilon = 0) \
&= [2 a(w, t) + 2 epsilon a(t, t)]_(epsilon = 0) \
&= 2 a(w, t) \
&= 2 thin (w, t) =^! 0 wide forall t in T_h
$ <eq:minimizer-perpencular>
This inner product being $0$ means that $w$ is perpendicular to $T_h$, or, $w in T_h^perp$.



#figure(
	cetz.canvas({
		import cetz.draw: line, circle, scale, content, arc

		scale(2.2)
		line((-2, 0), (2, 0), mark: (end: ">"))
		line((0, -2), (0, 2), mark: (end: ">"))
		content((-2, 2), anchor: "north-west", text(size: 18pt, $S_h$))
		let l = 1.2
		line((-l, -l), (l, l), stroke: 2pt)
		let o = 0.3
		line((l, l), (l + o, l + o), stroke: (dash: "dashed", thickness: 2pt))
		line((-l, -l), (-l - o, -l - o), stroke: (dash: "dashed", thickness: 2pt))
		content((l + o, l + o), anchor: "south-west", text(size: 18pt, $T_h$))
		
		let a = 1
		arc((a, a), start: -45deg, stop: -135deg, anchor: "origin", radius: 0.2)
		circle((a, a - 0.11), stroke: none, fill: black, radius: 0.5pt)
		line((0, 0), (a, a), mark: (end: ">"), stroke: (thickness: 2pt, paint: red), fill: red)
		let o = -0.5
		line((a, a), (a - o, a + o), mark: (end: ">"), stroke: (thickness: 2pt, paint: olive), fill: olive)
		line((0, 0), (a - o, a + o), mark: (end: ">"), stroke: (thickness: 2pt, paint: blue), fill: blue)

		content((a/2 + 0.2, a/2 + 0.2), $t = P u$, anchor: "south-east", padding: 0.1)
		content((a - o/2, a + o/2), $w = Q u = (I - P) u$, anchor: "south-west", padding: 0.1)
		content(((a - o) / 2, (a + o) / 2), $u$, anchor: "south-east", padding: 0.1)

	})
)

Since $G^"I"$ only changes the values at the black nodes, we can find a $t in T_h$ such that $u = w + t$. According to the orthogonal projection theorem, for a given vector $u$, vectors $w$ and $t$ are unique. Since $w in T_h^perp$ and $t in T_h$, we can view $t$ as the orthogonal projection of $u in S_h$ onto the subspace $T_h$. This is written as $P u = t$ where $P$ is the projection operator onto $T_h$. Similarly, we have $Q := I - P$ (where $I$ is the identity operator), such that
$
Q u = (I - P) u = I u - P u = u - t = w = G^"I" u.
$
Therefore, $Q = G^"I"$.

#hi[formal falsch, weil $G^"I"$ auf Vektoren arbeitet und $Q$ auf Funktionen (?)]

#line(length: 100%)

Now we want to show that $Q = G^"I"$ is self-adjoint, i.e.,
$
(Q u, v) = (u, Q v) wide forall u, v in S_h.
$

#hi[then be directly shown for $Q$ as a projection onto $T_h^perp$]

Let $u = u_1 + u_2$ where $u_1 in T_h$ and $u_2 in T_h^perp$. Similarly, $v = v_1 + v_2$ where $v_1 in T_h$ and $v_2 in T_h^perp$.
$
(P u, v) &= (P (u_1 + u_2), v) = (P u_1 + P u_2, v) = (u_1, v) = (u_1, v_1 + v_2) \
&= (u_1, v_1) + underbrace((u_1, v_2), =0 "since" u_1 perp v_2) = (u_1, v_1)
$
$
(u, P v) &= (u, P (v_1 + v_2)) = (u, P v_1 + P v_2) = (u, v_1) = (u_1 + u_2, v_1) \
&= (u_1, v_1) + underbrace((u_2, v_1), = 0 "since" u_2 perp v_1) = (u_1, v_1)
$

Therefore, $(P u, v) = (u, P v)$. The same is true for $Q$:
$
(Q u, v) = ((I - P) u, v) = (u - P u, v) = (u, v) - (P u, v) = (u, v) - (u, P v) \
= (u, v - P v) = (u, (I - P) v) = (u, Q v).
$ <eq:self-adjoint>

#line(length: 100%)

The same can be shown for $G^"II"$, only that here the relevant subspace is $hat(T)_h := {w in S_h | w(p_i) = 0 "for" p_i in Omega_h \\ Omega_H}$. Note that both $G^"I"$ and $G^"II"$ are idempotent. This means that if we apply the same orthogonal projection twice, the result will be the same as when we only apply it once ($G^"I" G^"I" = G^"I"$, or, $G^"I" G^"I" u = G^"I" u$). Assuming that $G^"II" u = u$,
$
||G^"I" u||^2 = (G^"I" u, G^"I" u) = (u, G^"I" G^"I" u) = (u, G^"I" u) \
= (G^"II" u, G^"I" u) = (u, G^"II" G^"I" u) <= ||u|| dot ||G^"II" G^"I" u||.
$
Hence,
$
(||G^"I" u||) / (||u||) <= (||G^"II" G^"I" u||) / (||G^"I" u||).
$
The same can be shown the other way around. Assuming $G^"I" u = u$,
$
(||G^"II" u||) / (||u||) <= (||G^"I" G^"II" u||) / (||G^"II" u||).
$
In other words, the norm reduction becomes less effective with each half-step. (Of course, the smaller $(||G^"I/II" u||) / (||u||)$ the better, since the solution to the homogeneous equation $laplace u = 0$ with a zero Dirichlet boundary condition is $u = 0$).

= Multigrid Algorithm

First iteration on level $q$ with $r = 2$ smoothings:

+ $u^(q, 0) <- G_q^"I" u^(q, 0)$ (special case because this is the first iteration)
+ $u^(q, 1) <- G_q^"II" u^(q, 0)$
+ $u^(q, 2) <- G_q^"I" u^(q, 1)$ (transition step)
+ approximate $v_1$ as the error $e in S_(q - 1)$ of the variational problem $J(u^(q, 2) + e) -> min$
	- on the coarsest level $q = 1$, compute $v_1$ exactly
	- otherwise, use a multigrid iteration on level $q - 1$, starting with $u^(q - 1, 0) = 0$
+ $u^(q, 3) <- u^(q, 2) + v_1$ (correction)
+ $u^(q, 5) <- G_q^"I" u^(q, 4)$ (mistake in the indices)
+ $u^(q, 6) <- G_q^"II" u^(q, 5)$
+ $u^(q, 7) <- G_q^"I" u^(q, 6)$
+ start the next iteration using $u^(q, 0) <- u^(q, 6)$ (even though we just computed $u^(q, 7)$)

#hi[
probably wrong: 6. should probably start with $G_q^"II"$
]

Next iterations omit step 1. and directly start with $u^(q, 1) = G_q^"II" u^(1, 0)$.

$
u^0 -->^(G^"I") u^1 -->^(G^"II") u^2 -->^(G^"I") u^3 -->^"coarse grid" u^4 -->^(G^"I") u^5 -->^(G^"II") u^6 -->^(G^"I") u^7
$

#line(length: 100%)

First iteration on level $q$ with $r = 3$ smoothings:

+ $u^(q, 0) <- G_q^"II" u^(q, 0)$ (special case because this is the first iteration)
+ $u^(q, 1) <- G_q^"I" u^(q, 0)$
+ $u^(q, 2) <- G_q^"II" u^(q, 1)$
+ $u^(q, 3) <- G_q^"I" u^(q, 2)$
+ coarse grid correction
+ $u^(q, 4) <- u^(q, 3) + v_1$
+ $u^(q, 6) <- G_q^"I" u^(q, 5)$ (mistake in the indices)
+ $u^(q, 7) <- G_q^"II" u^(q, 6)$
+ $u^(q, 8) <- G_q^"I" u^(q, 7)$
+ $u^(q, 9) <- G_q^"II" u^(q, 8)$
+ start the next iteration using $u^(q, 0) <- u^(q, 8)$

$
u^0 -->^(G^"II") u^1 -->^(G^"I") u^2 -->^(G^"II") u^3 -->^(G^"I") u^4 -->^"coarse grid" u^5 -->^(G^"I") u^6 -->^(G^"II") u^7 -->^(G^"I") u^8 -->^(G^"II") u^9
$

#line(length: 100%)

= Actual Proof

If we look at the sequence of Gauss-Seidel relaxations as an alternate method where $Q_V = G_q^"I"$ and $Q_W = G_q^"II"$,

#hi[why can we do that? How is $G_q^"I" = "id" - P_V$?]

then we notice that the last half-step is the least effective (since the norm reduction decreases from step to step):

$
(||u^(nu + 1)||) / (||u^(nu + 1/2)||) >= (||u^(nu + 1/2)||) / (||u^nu||).
$
This is independent of whether the last half-step is an application of $G_q^"I"$ or $G_q^"II"$ (since we can swap $Q_V$ and $Q_W$).


We also found that
$
(||G_q^"II" u^r||) / (||u^r||) <= (|u^r|) / (||u^r||) <= sqrt(rho).
$
This must be true even for the last application of $G_q^"II"$ in the sequence. All previous half-steps must have a norm reduction of $sqrt(q)$ or better (_better_ in this case means _smaller_ than $sqrt(q)$). Therefore,
$
(||u^r||) / (||u^0||) <= (sqrt(q))^r = q^(r/2)
$ <eq:6-2>

If $u' = u^r - u^(q - 1)$ then

#hi[Shouldn't it be $u' = u^r + u^(q - 1)$?]

#hi[proof]

$
||u'|| &<= sqrt(1/2 (1 - lambda^2)) dot ||u^r|| \
&<= rho^(r/2) sqrt((1 - lambda^2) / 2) dot ||u^0|| \
&= rho^(r/2) sqrt(((1 - lambda^2) / (lambda^2 + 1)) / (2 / (lambda^2 + 1))) dot ||u^0||
= rho^(r/2) sqrt(((lambda^2 + 1 - 2 lambda^2) / (lambda^2 + 1)) / ((2 lambda^2 + 2 - 2 lambda^2) / (lambda^2 + 1))) dot ||u^0||
= rho^(r/2) sqrt((1 - (2 lambda^2) / (lambda^2 + 1)) / (2 - (2 lambda^2) / (lambda^2 + 1))) dot ||u^0|| \
&= rho^(r/2) sqrt((1 - rho) / (2 - rho)) dot ||u^0||
$ <eq:6-3>

#infobox[
	#set math.equation(numbering: none)
	as a reminder:
	$
	rho := (2 lambda^2) / (lambda^2 + 1) wide lambda := (|v|) / (||v||)
	$
]

#line(length: 100%)

Now we look at the deviation of $v_1$ from $u^(q - 1)$. We find that
$
u^(r + 1) - u' = (u^r + v_1) - (u^r - u^(q - 1)) = v_1 + u^(q - 1).
$

In step 3 of the Multigrid algorithm, we require that
$
||u^(r + 1) - u'|| = ||v_1 + u^(q - 1)|| <^! delta ||u^(q - 1)||.
$ <eq:delta-estimation>

#hi[define $delta$ first]

We define
$
v_2 := 1/delta thin (u^(r + 1) - u') = 1/delta thin (v_1 + u^(q - 1)).
$

From @eq:delta-estimation, we find that
$
||u^(r + 1) - u'|| = delta ||v_2|| < delta ||u^(q - 1)|| \
==> ||v_2|| < ||u^(q - 1)||
$ <eq:v2vsuqm1>

Since $u'$ minimizes $J(v)$ in $u^r + S_(q - 1)$, we can show (similarly to @eq:minimizer-perpencular) that $u' perp u^(q - 1)$. Furthermore, we can apply the operator notation again: $P u^r = u^(q - 1) in S_(q - 1)$ is the projection of $u^r$ onto the subspace $S_(q - 1)$, and $Q u^r = (I - P) u^r = u^r - P u^r = u' in S_(q - 1)^perp$.

#figure(
	cetz.canvas({
		import cetz.draw: line, circle, scale, content, arc

		scale(2.2)
		line((-2, 0), (2, 0), mark: (end: ">"))
		line((0, -2), (0, 2), mark: (end: ">"))
		content((-2, 2), anchor: "north-west", text(size: 18pt, $S_q$))
		let l = 1.2
		line((-l, -l), (l, l), stroke: 2pt)
		let o = 0.3
		line((l, l), (l + o, l + o), stroke: (dash: "dashed", thickness: 2pt))
		line((-l, -l), (-l - o, -l - o), stroke: (dash: "dashed", thickness: 2pt))
		content((l + o, l + o), anchor: "south-west", text(size: 18pt, $S_(q-1)$))
		
		let a = 1
		arc((a, a), start: -45deg, stop: -135deg, anchor: "origin", radius: 0.2)
		circle((a, a - 0.11), stroke: none, fill: black, radius: 0.5pt)
		line((0, 0), (a, a), mark: (end: ">"), stroke: (thickness: 2pt, paint: red), fill: red)
		let o = -0.5
		line((a, a), (a - o, a + o), mark: (end: ">"), stroke: (thickness: 2pt, paint: olive), fill: olive)
		line((0, 0), (a - o, a + o), mark: (end: ">"), stroke: (thickness: 2pt, paint: blue), fill: blue)

		content((a/2 + 0.2, a/2 + 0.1), $u^(q - 1) = \ P u^r$, anchor: "south-east")
		content((a - o/2, a + o/2), $u' = Q u^r$, anchor: "south-west", padding: 0.1)
		content(((a - o) / 2, (a + o) / 2), $u^r$, anchor: "south-east")

	})
)

Since $v_2 in S_(q - 1)$, we find that $v_2 perp u'$, and therefore we can apply the Pythagorean theorem
$
||u' + v_2||^2 = ||u'||^2 + ||v_2||^2.
$
Because of @eq:v2vsuqm1, we have
$
||u' + v_2||^2 = ||u'||^2 + ||v_2||^2 <= ||u'||^2 + ||u^(q - 1)||^2 = ||u^r||^2.
$

#line(length: 100%)

In both the pre-smoothing and the post-smoothing, we have $r + 1$ alternating factors.
#figure(
	table(
		columns: 3,
		stroke: none,
		[*pre-smoothing*], [*coarse grid correction*], [*post-smoothing*],
		$underbrace(G^(r + 1) dot G^r dot ... dot G^2 dot G^1, =: tilde(G))$, $$, $underbrace(G^1 dot G^2 dot ... dot G^r dot G^(r + 1), = tilde(G)^*)$,
	)
)

We define
$
G^k := cases(G^"I" &wide& k "odd", G^"II" &wide& k "even")
$

Therefore,
$
u^r = tilde(G) u^0 wide "and" wide u^(2r + 2) = tilde(G)^* u^(r + 1).
$

Interestingly, $tilde(G)^*$ happens to be the adjoint of $tilde(G)$. This is because of the following reason:
$
(tilde(G) u, v) = (G^(r + 1) G^r G^(r - 1) ... G^1 u, v) &= (G^r G^(r - 1) ... G^1 u, G^(r + 1) v) \
&= (G^(r - 1) ... G^1 u, G^r G^(r + 1) v) \
&= ... \
&= (G^1 u, G^2 ... G^r G^(r + 1) v) \
&= (u, G^1 ... G^r G^(r + 1) v) = (u, tilde(G)^* v).
$
Note that $G^k$ are self-adjoint, as shown in @eq:self-adjoint.

For an arbitrary $hat(u) in S_q$, we find
$
(hat(u), u^(2 r + 2)) &= (hat(u), tilde(G)^* u^(r + 1)) \
&= (tilde(G) hat(u), u^(r + 1)) \
&= (tilde(G) hat(u), u' + delta v_2) \
&= (tilde(G) hat(u), u' - delta u' + delta u' + delta v_2) \
&= (tilde(G) hat(u), [1 - delta] u' + delta [u' + v_2]) \
&= [1 - delta]lr((tilde(G) hat(u), underbrace(u', Q dot Q u^r)), size: #40%) + delta (tilde(G) hat(u), u' + v_2) \
&= [1 - delta](Q tilde(G) hat(u), Q tilde(G) u^0) + delta (tilde(G) hat(u), u' + v_2) \
&<= [1 - delta](Q tilde(G) hat(u), Q tilde(G) u^0) + delta ||tilde(G) hat(u)|| dot ||u' + v_2|| \
&<= [1 - delta] ||Q tilde(G) hat(u)|| dot ||Q tilde(G) u^0|| + delta ||tilde(G) hat(u)|| dot ||u^r|| \
&= [1 - delta] ||Q tilde(G) hat(u)|| dot ||Q tilde(G) u^0|| + delta ||tilde(G) hat(u)|| dot ||tilde(G) u^0|| \
&<= sqrt((1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2) dot sqrt((1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2)
$ <eq:6-7>

The last transformation step can be explained as follows. If
$
a &:= sqrt(1 - delta) ||Q tilde(G) hat(u)||, &wide& c := sqrt(1 - delta) ||Q tilde(G) u^0||, \
b &:= sqrt(delta) ||tilde(G) hat(u)||, &wide& d := sqrt(delta) ||tilde(G) u^0||,
$
then
$
[1 - delta] ||Q tilde(G) hat(u)|| dot ||Q tilde(G) u^0|| + delta ||tilde(G) hat(u)|| dot ||tilde(G) u^0|| = a c + b d = vec(a, b) dot vec(c, d) \
<= lr(||vec(a, b)||) dot lr(||vec(c, d)||) = sqrt(a^2 + b^2) dot sqrt(c^2 + d^2) \
= sqrt((1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2) dot sqrt((1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2).
$

Let's focus on the term $(1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2$.

#infobox[
as a reminder:
#set math.equation(numbering: none)
$
||u^r|| &<= rho^(r / 2) ||u^0|| &wide& #ref(<eq:6-2>) \
||u'|| &<= rho^(r / 2) sqrt((1 - rho) / (2 - rho)) dot ||u^0|| &wide& #ref(<eq:6-3>)
$
]
From @eq:6-2 and @eq:6-3, we find that
$
(1 - delta) ||Q underbrace(tilde(G) u^0, u^r)||^2 + delta ||tilde(G) u^0||^2 &= (1 - delta) ||u'||^2 + delta ||u^r||^2 \
&<= (1 - delta) (rho^(r/2) sqrt((1 - rho) / (2 - rho)) ||u^0||)^2 + delta (rho^(r/2) ||u^0||)^2 \
&= (1 - delta) rho^r (1 - rho) / (2 - rho) ||u^0||^2 + delta rho^r ||u^0||^2 \
&= rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta) ||u^0||^2
$ <eq:6-8>

#infobox[
as a reminder:
$
rho = (2 lambda^2) / (lambda^2 + 1) wide "and" wide lambda =  (|v^r|) / (||v^r||)
$
]

Note that the starting vector $u^0$ in @eq:6-2, @eq:6-3 and @eq:6-8 can be set arbitrarily. However, starting from a different $u^0$, the resulting $u^r = tilde(G) u^0$ will be different as well, which leads to a different $rho$. Therefore, when we write @eq:6-8 for an arbitrary $u in S_q$ as a starting vector, we need to account for all possible values of $rho$.

#figure({
	let xs = lq.linspace(0, 4, num: 200)
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
})

Since $|v^r| <= ||v^r||$, we know that $0 <= lambda <= 1$. In this interval, $0 <= rho <= 1$. To preserve the inequality in @eq:6-8 for an arbitrary starting vector $u in S_q$, we need to include the upper bound of the term $rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)$ on this interval. Hence, we get
$
(1 - delta) ||Q tilde(G) u||^2 + delta ||tilde(G) u||^2 &<= underbrace(max_(0 <= rho <= 1) [rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)], =: delta_q) ||u||^2 \
&= delta_q ||u||^2
$

Therefore, we can extend equation @eq:6-7 as follows:

$
(hat(u), u^(2 r + 2)) &<= sqrt((1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2) dot sqrt((1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2) \
&<= sqrt(delta_q ||hat(u)||^2) dot sqrt(delta_q ||u^0||^2) \
&= delta_q dot ||hat(u)|| dot ||u^0|| wide forall hat(u) in S_q.
$

This inequality must be true for all $hat(u) in S_q$. In order to find an estimate for $||u^(2 r + 2)||$, we use the following trick:

$
(hat(u), u^(2 r + 2)) &<= delta_q ||hat(u)|| ||u^0|| \
=> ((hat(u), u^(2 r + 2))) / (||hat(u)||) &<= delta_q ||u^0|| \
=> (hat(u) / (||hat(u)||), u^(2 r + 2)) &<= delta_q ||u^0||.
$
The term $hat(u) / (||hat(u)||)$ is a normalized vector with length 1. Therefore, this inner product reaches its maximum, when $hat(u) = u^(2 r + 2)$. Since the inequality must hold for all $hat(u)$, it must hold for $hat(u) = u^(2 r + 2)$ as well:
$
max_(hat(u)) (hat(u) / (||hat(u)||), u^(2 r + 2)) &= (u^(2 r + 2) / (||u^(2 r + 2)||), u^(2 r + 2)) = ((u^(2 r + 2), u^(2 r + 2))) / (||u^(2 r + 2)||)  = (lr(||u^(2 r + 2)||)^2) / (||u^(2 r + 2)||) \
&= ||u^(2 r + 2)|| <= delta_q ||u^0||.
$

Now we have an estimate of the error norm at the end of a Multigrid cycle compared to its norm at the beginning of the cycle.

#line(length: 100%)

= Convergence of the Two-Grid Method

When we only have two different levels, we observe the following: On the coarser level, we solve for the error in step 3 of the Multigrid algorithm exactly. Therefore,
$
0 = ||v_1 - u^(q - 1)|| < delta ||u^(q - 1)|| wide ==> wide delta = 0.
$

In this case,
$
delta_"TG" := delta_q = max_(0 <= rho <= 1) [rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)] = max_(0 <= rho <= 1) underbrace([rho^r (1 - rho) / (2 - rho)], delta_"TG" (rho)).
$

#let maxima = range(5).map(r => {
	if r == 0 {
		return (0, 0.5)
	}
	let a = (3 * r + 1) / (2 * r)
	let b = a - calc.sqrt(a * a - 2)
	let c = calc.pow(b, r) * (1 - b) / (2 - b)
	return (b, c)
})

#figure(
	lq.diagram(
		..range(5).map(r => lq.plot(
			lq.linspace(1e-8, 1),
			x => calc.pow(x, r) * (1 - x) / (2 - x),
			mark: none,
			label: $r = #r$,
		)),
		lq.scatter(
			maxima.map(m => m.at(0)),
			maxima.map(m => m.at(1)),
			color: black,
		),
		legend: (
			position: (100% + 0.5em, 0%),
		),
		xaxis: (
			label: $rho$,
		),
		yaxis: (
			label: $delta_"TG" (rho)$
		),
	)
)

The maxima of $rho_"TG" (rho)$ can be found analytically by setting its derivative to zero. However, the resulting term for the maximum is rather complicated. Braess found a good upper bound:
$
rho_"TG" < 1/((r + 1) e).
$

= Convergence of the Multigrid Method

#hi[
Convergence rate on level $q - 1$ must be the same as on level $q$ when $q > 1$.
]

In case $delta = delta_("MG", mu)^mu$, we find
$
delta_("MG", mu) &= max_(0 <= rho <= 1) rho^r ((1 - delta_("MG", mu)^mu) (1 - rho) / (2 - rho) + delta_("MG", mu)^mu) \
&= max_(0 <= rho <= 1) rho^r ((1 - rho) / (2 - rho) + delta_("MG", mu)^mu (1 - (1 - rho) / (2 - rho))) \
&= max_(0 <= rho <= 1) rho^r ((1 - rho) / (2 - rho) + delta_("MG", mu)^mu / (2 - rho)) \
&= max_(0 <= rho <= 1) rho^r ((1 + delta_("MG", mu)^mu - rho) / (2 - rho)) \
$

This equation can be solved numerically.