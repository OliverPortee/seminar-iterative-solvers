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

#let Gbox(content, color) = box(fill: color.transparentize(70%), inset: (top: 7pt, bottom: 4pt, left: 3pt, right: 2pt), radius: 5pt, stroke: color, baseline: 15%, content)
#let G1(content) = Gbox(content, fine-color)
#let G2(content) = Gbox(content, coarse-color)

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
	caption: [Each shaded area corresponds to one of the basis functions of $T_h$. The supports of the respective basis functions are disjoint (they do not overlap each other).]
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
sum_nu C_nu = sum_vec(i\, j in Omega_(2h), d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2
$

We can also show that
$
sum_nu D_nu >= sum_vec(i\, j in Omega_(2h) without Omega_H, d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2
$

Together with $C$, we get
$
sum_nu C_nu + D_nu &>= sum_vec(i\, j in Omega_(2h), d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2 + sum_vec(i\, j "red", d(i, j) = 2 h, delim: #none) 1/4 (v_i - v_j)^2 \
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
Since $u$ and $v$ have the same values on the coarse points, we know that $|u| = |v|$.
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
	d(p1, blue)
	d(p2, yellow)
	d(p3, blue)
	d(p4, yellow)
	d(p5, blue)
	d(p6, yellow)
	d(p7, blue)
	d(p8, yellow)
	d(p0, blue)
	point-labels()
	content(p0, $0$, anchor: "north-west", padding: 0.05)
	
	translate(x: 6)
	
	coordinate-system()
	line(p1, p3, p5, p7, close: true)
	line(p3, p7)
	line(p1, p5)
	line(p4, p8)
	line(p2, p6)
	d(p1, yellow)
	d(p2, blue)
	d(p3, yellow)
	d(p4, blue)
	d(p5, yellow)
	d(p6, blue)
	d(p7, yellow)
	d(p8, blue)
	d(p0, yellow)
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
(G_h^"I" u)_i &= cases(u_i &wide& p_i in Omega_H wide ("yellow") &wide& , 1/4 sum_(j in N(i)) u_j &wide& p_i in Omega_h without Omega_H wide ("blue") 
) \
(G_h^"II" u)_i &= cases(1/4 sum_(j in N(i)) u_j &wide& p_i in Omega_H wide("yellow") &wide& , u_i &wide& p_i in Omega_h without Omega_H wide ("blue") )
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
		d(p1, yellow)
		d(p2, blue)
		d(p3, yellow)
		d(p4, blue)
		d(p5, yellow)
		d(p6, blue)
		d(p7, yellow)
		d(p8, blue)
		d(p0, yellow)
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

