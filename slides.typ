
#import "@preview/fletcher:0.5.8" as fl
#import "@preview/touying:0.7.4": *
#import themes.simple: *

#import "lib.typ": *

#show: simple-theme

#let hi(content) = box(width: 100%, inset: 0.8em, stroke: orange, content)

#title-slide[
#v(1fr)
= MG Konvergenzrate
Teil 2

#v(1fr)
Oliver Portee
]

= Gauß-Seidel Relaxation

#let Gbox(content, color) = box(fill: color.transparentize(70%), inset: (top: 7pt, bottom: 4pt, left: 3pt, right: 2pt), radius: 5pt, stroke: color, baseline: 15%, content)
#let G1(content) = Gbox(content, blue)
#let G2(content) = Gbox(content, olive)

#figure(
	cetz.canvas({
		import cetz.draw: line, circle, scale

		let grid-point(x, y) = {
			let fill = if even(x + y) {
				olive
			} else {
				blue
			}
			cetz.draw.circle((x, y), fill: fill, stroke: none, radius: 4pt)
		}
		
		scale(2)

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
	})
)


#G1($G^"I"$) glättet die blauen Punkte #h(1fr) #G2($G^"II"$) glättet die grünen Punkte


= Multigrid-Algorithmus


// #fl.diagram(
// 	spacing: (1.5em, 2em),
// 	fl.node((0, 0), $u^"start"$),
// 	fl.edge($G_(r + 1)$, "-|>", bend: 40deg),
// 	fl.node((1, 0), $u^1$),
// 	fl.edge($G_r$, "-|>", bend: 40deg),
// 	fl.node((2, 0), $...$),
// 	fl.edge($G_1$, "-|>", bend: 40deg),
// 	fl.node((3, 0), $u^r$),
// 	fl.node((4, 0), [coarse grid \ correction]),
// 	fl.node((5, 0), $u^(r + 1)$),
// 	fl.edge($G_1$, "-|>", bend: 40deg),
// 	fl.node((6, 0), $u^(r + 2)$),
// 	fl.edge($G_2$, "-|>", bend: 40deg),
// 	fl.node((7, 0), $...$),
// 	fl.edge($G_(r + 1)$, "-|>", bend: 40deg),
// 	fl.node((8, 0), $u^(2 r + 2)$),
// )

#let fine(content) = $markhl(content)$

#figure(
	cetz.canvas({
		import cetz.draw: line, content, scale

		scale(2.2)
		let r = 4
		for x in range(r + 1) {
			if x > 0 {
				line((x - 0.70, 0), (x - 0.3, 0), mark: (end: ">"))
			}
			if calc.even(x) {
				content((x, 0), G1($G^#(r + 1 - x)$))
			} else {
				content((x, 0), G2($G^#(r + 1 - x)$))
			}
		}
		line((r + 0.1, -0.3), (r + 0.5, -1.0), mark: (end: ">"))
		content((r + 1, -1.5), align(center)[coarse grid \ correction])
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
	})
)

#v(2em)


$
G^k := cases(#G1($G^"I"$) &wide& k "odd", #G2($G^"II"$) &wide& k "even")
$

== Was haben wir bereits gezeigt?

$
	|u| <= sqrt((2 lambda^2) / (1 + lambda^2)) ||u|| quad "wobei" quad lambda = (|v|) / (||v||)
$

$
	||G_h^"II" u||^2 <= |u| quad "falls" quad u = G_h^"I" u
$

Wo wollen wir hin?

$
	||u^(2 r + 2)|| <= delta_q ||u^0||
$

== Herangehensweise

+ Konvergenzrate nach $r + 1$ Gauß-Seidel-Halbschritten
+ Konvergenzrate nach Grobgitterkorrektur
+ ...

= Konvergenzrate beim Pre-Smoothing


== Aufspaltung des Funktionenraums
$
S_h = S_H plus.o T_h
$

$S_h$ finite Elemente Raum mit einem Freiheitsgrad an jedem Gitterpunkt

$S_H$ finite Elemente Raum mit einem Freiheitsgrad an jedem Grobgitterpunkt

$T_h = {w in S_h | w(p_i) = 0 "falls" p_i #G2("Grobgitterpunkt")}$

#hi[Visualization of $T_h$]

== Gauß-Seidel-Halbschritte als Projektionen

$
(G_h^"I" u)_i := cases(1/4 sum_(j in N(i)) u_j &wide& p_i in Omega_h \\ Omega_H, u_i &wide& p_i in Omega_H)
$

#figure(
	cetz.canvas({
		import cetz.draw: circle, line, content, scale, floating

		let grid-point(coords) = circle(coords, radius: 2pt, stroke: none, fill: if even(coords.at(0) + coords.at(1)) { blue } else { olive } )

		let p0 = ( 0,  0)
		let p1 = ( 0, -1)
		let p2 = ( 1,  0)
		let p3 = ( 0,  1)
		let p4 = (-1,  0)
		let p5 = ( 1, -1)
		let p6 = ( 1,  1)
		let p7 = (-1,  1)
		let p8 = (-1, -1)

		scale(3)

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

		floating({
			content(p0, $0$, anchor: "north-west", padding: 0.3)
			content(p1, $1$, anchor: "north", padding: 0.3)
			content(p2, $2$, anchor: "west", padding: 0.3)
			content(p3, $3$, anchor: "south", padding: 0.3)
			content(p4, $4$, anchor: "east", padding: 0.3)

			content((0, -1.6), $u_0 <--^(G_h^"I") 1/4 (u_1 + u_2 + u_3 + u_4)$)
		})
	})
)


== Gauß-Seidel-Halbschritte als Projektionen

$G_h^"I"$ minimiert $J(v)$ über $u + T_h$, denn:
$
J(v) :&= a(v, v) - cancel(2 (f, v)_0) \
&= ||v||^2 \
&= sum_vec(i\, j, d(i, j) = h, delim: #none) (v_i - v_j)^2
$

== Gauß-Seidel-Halbschritte als Projektionen

$G_h^"I"$ minimiert $J(v)$ über $u + T_h$, denn:
$
v_0 = 1/4 (v_1 + v_2 + v_3 + v_4) \
"minimiert" \
(v_1 - v_0)^2 + (v_2 - v_0)^2 + (v_3 - v_0)^2 + (v_4 - v_0)^2
$

== Gauß-Seidel-Halbschritte als Projektionen

$
v := G_h^"I" u \
J(v + epsilon w) >= J(u) wide forall w in T_h, epsilon in RR
$

alternativ:
$
lr(dv(, epsilon) J(v + epsilon w)|)_(epsilon = 0) = 0 wide forall w in T_h
$

== Gauß-Seidel-Halbschritte als Projektionen

$
dv(, epsilon) J(v + epsilon w) &= dv(, epsilon) a(v + epsilon w, v + epsilon w) \
&= dv(, epsilon) [a(v, v) + 2 epsilon a(v, w) + epsilon^2 a(w, w)] \
&= 2a(v, w) + 2 epsilon a(w, w)
$
für $epsilon = 0$:
$
a(v, w) =^! 0 wide forall w in T_h
$

== Gauß-Seidel-Halbschritte als Projektionen

#figure(
	cetz.canvas({
		import cetz.draw: line, content, scale, floating, circle
		scale(1.8)
		line((-1, -1), (5, 5), stroke: 2pt)
		floating({
			content((-4, 5), $S_h = S_H plus.o T_h$, anchor: "west")
			content((5, 5), $T_h$, anchor: "west", padding: 0.3)
			content((3, 0), $u in S_h$, anchor: "north", padding: 0.2)
			content((4.5, 1.5), $v = G_h^"I" u in T_h^perp$, anchor: "south-west", padding: 0.2)
			content((1.5, 1.5), $w in T_h$, anchor: "south-east", padding: 0.2)
		})
		circle((0, 0), radius: 3pt, stroke: none, fill: black)
		content((0, 0), $va(0)$, anchor: "south-east", padding: 0.2)
		line((0, 0), (6, 0), mark: (end: ">"), stroke: 5pt)
		line((3, 3), (6, 0), mark: (end: ">"), stroke: 5pt)
		line((0, 0), (3, 3), mark: (end: ">"), stroke: 5pt)
	})
)

== Gauß-Seidel-Halbschritte als Projektionen

- sei $P u$ die Orthogonalprojektion von $u$ auf $T_h$
#pause
- dann ist $w = P u$ und $v = u - P u = underbrace((I - P), G_h^"I") u$
#pause
- $G_h^"I"$ ist eine Orthogonalprojektion $T_h^perp$
#pause
- übertragbar auf $G_h^"II"$, bloß dass dann \
	$T_h = {w in S_h | w(p_i) = 0 "falls" p_i in Omega_h \\ Omega_H}$
#pause
- Orthogonalprojektionen sind selbstadjungiert

== Adjungierte Operatoren

Operator $Q^*$ ist der adjungierte Operator zu $Q$, wenn gilt
$
(Q u, v) = (u, Q^* v).
$

$Q$ ist selbstadjungiert, wenn gilt
$
Q^* = Q,
$
also
$
(Q u, v) = (u, Q v).
$

== Orthogonalprojektionen sind selbstadjungiert

Sei $u = u_1 + u_2 in U$ und $v = v_1 + v_2 in V$
