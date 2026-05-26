#import "lib.typ": *

#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")

#title[Convergence Rate of a Multigrid Method with Gauss-Seidel Relaxation for the Poisson Equation]

= Grids

#figure(
	cetz.canvas({

		import cetz.draw: *
		scale(0.7)
		scale(y: -100%)
		h-grid()
		content((W / 2, H + 1), $Omega_h$)
		translate(x: 8)
		H-grid()
		content((W / 2, H + 1), $Omega_H$)
		translate(x: 8)
		h2-grid()
		content((W / 2, H + 1), $Omega_(2 h)$)
	})
)

= Decomposition

#figure(
	cetz.canvas({
		import cetz.draw: *
		scale(0.7)
		scale(y: -100%)
		h-grid()
		content((W / 2, H + 1), $u in S_h$)
		content((7, H / 2), $=$)
		translate(x: 8)
		H-grid()
		content((W / 2, H + 1), $v in S_H$)
		content((7, H / 2), $plus.o$)
		translate(x: 8)
		h-grid(hide-green-red: true)
		content((W / 2, H + 1), $w in T_h$)
	})
)

$
S_h &= S_H plus.o T_h wide u &= v plus w
$

= Lemma 3.1

$
integral_I v_xi w_xi + v_eta w_eta = integral_I v_xi w_xi + integral_I v_eta w_eta
$

$
integral_I v_xi w_xi = integral_(I_"links") v_xi w_xi + integral_(I_"rechts") v_xi w_xi = 0
$
da $lr(w_xi|)_(I_"links") = - lr(w_xi|)_(I_"rechts")$ und $v_xi = "const"$.

$
integral_I v_xi w_xi + v_eta w_eta &= integral_I v_eta w_eta \
&= integral_I 1/2 v_eta^2 - 1/2 v_eta^2 + 1/2 w_eta^2 - 1/2 w_eta^2 + v_eta w_eta \
&= integral_I 1/2 (-v_eta^2 + 2 v_eta w_eta - w_eta^2) + 1/2 v_eta^2 + 1/2 w_eta^2 \
&= integral_I 1/2 v_eta^2 + 1/2 w_eta^2 - 1/2 (v_eta^2 - 2 v_eta w_eta + w_eta^2) \
&= 1/2 integral_I v_eta^2 + 1/2 integral_I w_eta^2 - 1/2 integral_I (v_eta + w_eta)^2 \
&= 1/2 integral_I v_eta^2 + 1/2 integral_I v_xi^2 - 1/2 integral_I v_xi^2 + 1/4 integral_I (w_eta^2 + w_eta^2) - 1/2 integral_I (v_eta + w_eta)^2 \
&= underbrace(1/2 integral_I (v_eta^2 + v_xi^2), A_I) + underbrace(1/4 integral_I (w_eta^2 + w_xi^2), B_I) - underbrace(1/2 integral_I v_xi^2, C_I) - underbrace(1/2 integral_I (v_eta + w_eta)^2, D_I) \
&= A_I + B_I - C_I - D_I
$

#highlight[TODO]

#figure(
	cetz.canvas({
		import cetz.draw: line, content, scale, circle

		scale(2)

		let p1 = (0, -1)
		let p2 = (1, 0)
		let p3 = (0, 1)
		let p4 = (-1, 0)
		let p5 = (0, 0)
		
		line((p4.at(0) - 0.5, 0), (p2.at(0) + 0.5, 0), mark: (end: ">"), name: "xi-axis")
		line((0, p1.at(1) - 0.3), (0, p3.at(1) + 0.4), mark: (end: ">"), name: "eta-axis")
		content("xi-axis.end", $xi$, anchor: "north-west", padding: 0.1)
		content("eta-axis.end", $eta$, anchor: "south-east", padding: 0.1)

		line(p1, p2, p3, p4, close: true)

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
	})
)

$
a(v, w) <= underbrace(1/2 ||v||^2, A) + underbrace(1/4 ||w||^2, B) - underbrace(1/2 |v|^2, C + D) = 1/2 (||v||^2 - |v|^2) + 1/4 ||w||^2
$

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