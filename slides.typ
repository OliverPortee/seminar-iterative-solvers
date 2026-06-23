
#import "@preview/fletcher:0.5.8" as fl
#import "@preview/touying:0.7.4": *
#import "@preview/cetz:0.5.2"
#import "@preview/physica:0.9.8": *
#import "@preview/lilaq:0.6.0" as lq

#import calc: even, odd
#import themes.simple: *


#show: simple-theme

#let hi(content) = box(width: 100%, inset: 0.8em, stroke: orange, content)

#title-slide[
#v(1fr)
= MG Konvergenzrate

Über die Kunst, Probleme so lange zu projizieren, bis sie konvergieren

Teil 2

#v(1fr)
Oliver Portee
]

= Gauß-Seidel Relaxation

#let fine-color = lq.color.map.petroff10.at(0)
#let coarse-color = lq.color.map.petroff10.at(1)

#let Gbox(content, color) = box(fill: color.transparentize(70%), inset: (top: 7pt, bottom: 4pt, left: 3pt, right: 2pt), radius: 5pt, stroke: color, baseline: 15%, content)
#let G1(content) = Gbox(content, fine-color)
#let G2(content) = Gbox(content, coarse-color)

#let draw-grid(show-th: false) = {
	import cetz.draw: line, circle

	let W = 6
	let H = 4

	let grid-point(x, y) = {
		let fill = if even(x + y) {
			coarse-color
		} else {
			fine-color
		}
		cetz.draw.circle((x, y), fill: fill, stroke: none, radius: 4pt)
	}

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

#figure(
	cetz.canvas({
		cetz.draw.scale(2)
		draw-grid()
	})
)


#G1($G^"I"$) glättet die blauen Punkte #h(1fr) #G2($G^"II"$) glättet die gelben Punkte


= Multigrid-Algorithmus


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
G^k := cases(#G1($G^"I"$) &wide& k "ungerade", #G2($G^"II"$) &wide& k "gerade")
$


==  Wo wollen wir hin?

$
	||u^(2 r + 2)|| <= delta_q ||u^0||
$

mit
$
delta_q < 1
$

== Herangehensweise

+ Konvergenzrate nach $r$ Gauß-Seidel-Halbschritten
+ Konvergenzrate nach Grobgitterkorrektur
+ Konvergenzrate nach einer vollständigen MG Iteration

= Konvergenzrate beim Presmoothing


== Aufspaltung des Funktionenraums
$
S_h = S_H plus.o T_h
$

$S_h$ finite Elemente Raum mit einem Freiheitsgrad an jedem Gitterpunkt

$S_H$ finite Elemente Raum mit einem Freiheitsgrad an jedem Grobgitterpunkt

$T_h = {w in S_h | w(p_i) = 0 "falls" p_i #G2("Grobgitterpunkt")}$


== Basisfunktionen von $T_h$

#slide(repeat: 2, self => [
	#figure(
		cetz.canvas({
			cetz.draw.scale(2)
			draw-grid(show-th: self.subslide > 1)
		})
	)

	$
	T_h = {w in S_h | w(p_i) = 0 "falls" p_i #G2("Grobgitterpunkt")}
	$
])

== Gauß-Seidel-Halbschritte als Projektionen

$
(G_h^"I" u)_i := cases(1/4 sum_(j in N(i)) u_j &wide& p_i in Omega_h \\ Omega_H, u_i &wide& p_i in Omega_H)
$

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
J(v + w) >= J(v) wide forall w in T_h
$

alternativ:
$
lr(dv(, epsilon) J(v + epsilon w)|)_(epsilon = 0) = 0 wide forall w in T_h, epsilon in RR
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
			content((1.5, 1.5), $P_(T_h) u = w in T_h$, anchor: "south-east", padding: 0.2)
		})
		circle((0, 0), radius: 3pt, stroke: none, fill: black)
		content((0, 0), $va(0)$, anchor: "south-east", padding: 0.2)
		line((0, 0), (6, 0), mark: (end: ">"), stroke: 5pt)
		line((3, 3), (6, 0), mark: (end: ">"), stroke: 5pt)
		line((0, 0), (3, 3), mark: (end: ">"), stroke: 5pt)
	})
)

#highlight[
TODO: genau begründen, warum $G^"I" u in T^perp$

1. $v = G^"I" u$
2. $v perp T_h ==> v in T_h^perp$
3. $T_h plus.o T_h^perp = S_h$
4. $u in S_h ==> exists w in T_h$ s. th. $u = w + v$
5. laut ... Satz ist diese Zerlegung eindeutig
]

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

#highlight[TODO: weniger Stichpunkte]

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

Sei $x = x_1 + x_2$ und $y = y_1 + y_2$, wobei $x_1, y_1 in T_h$ und $x_2, y_2 in T_h^perp$.

$
(P x, y) = (P x_1 + cancel(P x_2), y) = (x_1, y) = (x_1, y_1 + y_2) \
= (x_1, y_1) + cancel((x_1, y_2)) = (x_1, y_1) = (x_1, y_1) + cancel((x_2, y_1)) \
= (x_1 + x_2, y_1) = (x, y_1) = (x, P y_1 + cancel(P y_2)) = (x, P y)
$

== Projektionen sind idempotent

$
P^2 x = P P x = P x wide ==> wide P^2 = P
$

== Gauß-Seidel-Halbschritte werden ineffektiver

Falls $G_h^"II" u = u$, dann

$
||G_h^"I" u||^2 = (G_h^"I" u, G_h^"I" u) = (u, G_h^"I" G_h^"I" u) = (u, G_h^"I" u) \
= (G_h^"II" u, G_h^"I" u) = (u, G_h^"II" G_h^"I" u) <= ||u|| dot ||G_h^"II" G_h^"I" u||.
$

$
(||G_h^"I" u||) / (||u||) <= (||G_h^"II" G_h^"I" u||) / (||G_h^"I" u||)
quad "analog:" quad
(||G_h^"II" u||) / (||u||) <= (||G_h^"I" G_h^"II" u||) / (||G_h^"II" u||)
$

== Was haben wir bereits in Teil 1 gezeigt?


$
	||G_h^"II" u|| <= |u| quad "falls" quad u = G_h^"I" u
$

$
	|u| <= sqrt((2 lambda^2) / (1 + lambda^2)) ||u|| quad "wobei" quad lambda = (|v|) / (||v||)
$

#pause

$
==> ||G_h^"II" u|| <= |u| <= sqrt((2 lambda^2) / (1 + lambda^2)) ||u|| \
==> (||G_h^"II" u||) / (||u||) <= sqrt((2 lambda^2) / (1 + lambda^2)) =: sqrt(rho)
$

== Der letzte Halbschritt ist am ineffektivsten



#figure(
	{
		show lq.selector(lq.legend): set grid(row-gutter: 1em)
		let n = 10
		let rhoval = 1.03
		lq.diagram(
			height: 84%,
			width: 87%,
			legend: (
				inset: 13pt,
			),
			xaxis: (
				label: [Iteration],
				subticks: none,
			),
			yaxis: (
				subticks: none,
				ticks: none,
				extra-ticks: (lq.tick(rhoval, label: $sqrt(rho)$),),
				lim: (-0.05, auto),
			),
			lq.scatter(
				range(n),
				range(n).map(x => 2 * calc.exp(-x)),
				size: 23pt,
				label: $||u^k||$,
			),
			lq.scatter(
				range(n),
				range(n).map(x => 1 - 0.8 * calc.exp(-x)),
				size: 23pt,
				label: $(||u^(k + 1)||) / (||u^k||)$,
			),
			lq.line((0%, rhoval), (100%, rhoval)),
		)
	}
)

== Bei $r$ Halbschritten

$
||u^r|| <= (sqrt(rho))^r ||u^0|| = rho^(r/2) ||u^0||
$

$rho^(r / 2)$ ... Konvergenzrate des Presmoothers #emoji.checkmark.box

#highlight[TODO: Abbildung mit Smoothing-Steps]


= Konvergenzrate der Grobgitterkorrektur


Angenommen, we kennen schon die Konvergenzrate $delta_(q - 1)$:
$
delta := cases(0 &wide& "gröbstes Level", (delta_(q - 1))^mu &wide& "feinere Level")
$

$
mu := cases(1 &wide& "V-cycle", 2 &wide& "W-cycle")
$

$
||v_1 - u^(q - 1)|| < delta ||u^(q - 1)||
$


== Exakte Grobgitterkorrektur

$
u' = u^r - u^(q - 1)
$

Exakte Lösung:
$
u^(q - 1) := limits("argmin")_(u in S_H) space J(u^r - u)
$

Exakte Grobgitterkorrektur:
$
u' = u^r - u^(q - 1)
$

== Exakte Grobgitterkorrektur

Aus Teil 1:

$
a(v, w) <= sqrt(1/2 (1 - lambda^2)) ||v|| ||w|| wide forall v in S_H, w in T_h
$

zu zeigen:

$
||u'|| <= sqrt(1/2 (1 - lambda^2)) ||u^r||
$

#let picture = cetz.canvas({
		import cetz.draw: *

		let right-angle(pos, start: 0deg) = {
			arc(pos, anchor: "origin", start: start, stop: start + 90deg)
			let r = 0.57
			let dx = r * calc.cos(start + 45deg)
			let dy = r * calc.sin(start + 45deg)
			let label-pos = (pos.at(0) + dx, pos.at(1) + dy)
			circle(label-pos, stroke: none, fill: black, radius: 1pt)
		}
		
		scale(1.5)
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
		floating(content((15/2, 6), $T_h$, anchor: "south-east", padding: 0.1))

		arc((0, 0), anchor: "origin", start: 90deg, stop: 128deg)
		content((-0.21, 0.7), $alpha$)
		arc((0, 0), anchor: "origin", start: 0deg, stop: 38deg)
		content((0.7, 0.29), $alpha$)

		right-angle((0, 5), start: 180deg)
	})

#align(center + horizon,
	picture
)
$
u' = u^r - u^(q - 1) = u^r - P_(S_H) u^r
$


== Exakte Grobgitterkorrektur

#scale(68%, picture, reflow: true)

#place(top + left, dx: 446pt, dy: 34pt, text(20pt)[
	$
		v in S_H, w in T_h \
		cos(alpha) = a(v, w) / (||v|| ||w||) \
		a(v, w) <= sqrt(1/2 (1 - lambda^2)) ||v|| ||w|| \
		a(v, w) / (||v|| ||w||) = cos(alpha) <= sqrt(1/2 (1 - lambda^2)) \
	$
])

$
	||u'|| = ||u^r|| cos(alpha) <= sqrt(1/2 (1 - lambda^2)) ||u^r||
$

== Presmoothing + Exakte Grobgitterkorrektur

$
||u'|| <= sqrt(1/2 (1 - lambda^2)) ||u^r||  wide "und" wide ||u^r|| <= rho^(r / 2) ||u^0|| \
$
$
==> ||u'|| &<= rho^(r / 2) sqrt(1/2 (1 - lambda^2)) ||u^0|| \
&= rho^(r/2) sqrt((1 - rho) / (2 - rho)) ||u^0||
$
mit $rho := (2 lambda^2) / (lambda^2 + 1)$

#highlight[TODO: Abbildung mit Smoothing-Steps und Grobgitterkorrektur]

#highlight[possibly focus more on the first part]

== Presmoothing + Exakte Grobgitterkorrektur

#{
let xs = lq.linspace(0, 1, num: 250)
let r = 3
show lq.selector(lq.legend): set grid(row-gutter: 1em)
figure(
	lq.diagram(
		width: 100%,
		height: 84%,
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
	)
)
}

#place(
	top + left,
	dx: 79pt,
	dy: 85pt,
	$
		(r = 3) \
	$
)


== Exakte Grobgitterkorrektur ist Orthogonalprojektion

#figure(
	cetz.canvas({
		import cetz.draw: *

		scale(1.5)
		
		line((-3, 0), (8, 0))

		line((0, 0), (6, 6), stroke: 2pt, mark: (end: ">"))
		line((6, 0), (6, 6), stroke: 2pt, mark: (end: ">"))
		line((0, 0), (6, 0), stroke: 2pt, mark: (end: ">"))

		floating({
			content((3, 3), $u^r$, anchor: "south-east")
			content((3, 0), $u^(q - 1)$, anchor: "north", padding: 0.3)
			content((0, 0), $0$, anchor: "north", padding: 0.2)
			content((6, 3), $u' &= u^r - u^(q - 1) \ &= underbrace((I - P_(S_H)), =: Q) u^r$, anchor: "west", padding: 0.2)
			content((8, 0), $S_H$, anchor: "west", padding: 0.2)
		})
	})
)

== Nicht-Exakte Grobgitterkorrektur

$
||u^(r + 1) - u'|| = ||(u^r + v_1) - (u^r + u^(q - 1))|| = ||v_1 - u^(q - 1)|| \
$
laut Annahme:
$
||underbrace(v_1 - u^(q - 1), =: delta v_2)|| <= delta ||u^(q - 1)|| \
$
$
v_1 - u^(q - 1) = u^(r + 1) - u' =: delta v_2 wide "mit" wide ||v_2|| <= ||u^(q - 1)||
$


== Nicht-Exakte Grobgitterkorrektur

#figure(
	cetz.canvas({
		import cetz.draw: *

		scale(1.5)
		
		line((-3, 0), (8, 0))

		line((0, 0), (6, 6), stroke: 2pt, mark: (end: ">"))
		line((6, 0), (6, 6), stroke: 2pt, mark: (end: ">"))
		line((0, 0), (6, 0), stroke: 2pt, mark: (end: ">"))
		line((0, 0), (1.5, 0), stroke: 3pt, mark: (end: ">"))

		line((6, 6), (7.5, 6), stroke: (dash: "dashed"))
		line((6, 0), (7.5, 6), stroke: 2pt, mark: (end: ">"))
		line((1.5, 0), (7.5, 6), stroke: (dash: "dashed"))

		floating({
			content((3, 3), $u^r$, anchor: "south-east")
			content((6, 0), $u^(q - 1)$, anchor: "north", padding: 0.3)
			content((0, 0), $0$, anchor: "north", padding: 0.2)
			content((6, 3), $u'$, anchor: "east", padding: 0.1)
			content((6.75, 3), $u^(r + 1)$, anchor: "west", padding: 0.2)
			content((8, 0), $S_H$, anchor: "west", padding: 0.2)
			content((1.5, 0), $delta v_2$, anchor: "north", padding: 0.2)
		})
	})
)

#highlight[TODO: $G$ einführen]

== Fast Geschafft

$
(hat(u), u^(2 r + 2)) = (hat(u), tilde(G)^* u^(r + 1)) = (tilde(G) hat(u), u' + delta v_2) \
= (tilde(G) hat(u), u' -delta u' + delta u' + delta v_2) = (tilde(G) hat(u), (1 - delta) u' + delta (u' + v_2)) \
= (1 - delta) (tilde(G) hat(u), u') + delta (tilde(G) hat(u), u' + v_2) \
<= (1 - delta) (tilde(G) hat(u), Q u^r) + delta ||tilde(G) hat(u)|| ||u' + v_2|| \
<= (1 - delta) (tilde(G) hat(u), Q^2 u^r) + delta ||tilde(G) hat(u)|| ||u^r|| \
<= (1 - delta) ||Q tilde(G) hat(u)|| ||Q tilde(G) u^0|| + delta ||tilde(G) hat(u)|| ||tilde(G) u^0||
$

$
(hat(u), u^(2 r + 2)) &<= (1 - delta) ||Q tilde(G) hat(u)|| ||Q tilde(G) u^0|| + delta ||tilde(G) hat(u)|| ||tilde(G) u^0|| \
&= a c + b d = vec(a, b) dot vec(c, d) \
&<= lr(||vec(a, b)||) dot lr(||vec(c, d)||) = sqrt(a^2 + b^2) dot sqrt(c^2 + d^2) \
&= sqrt((1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2) dot sqrt((1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2)
$

$
a &:= sqrt(1 - delta) ||Q tilde(G) hat(u)|| &wide c &:= sqrt(1 - delta) ||Q tilde(G) u^0|| \
b &:= sqrt(delta) ||tilde(G) hat(u)|| &wide d &:= sqrt(delta) ||tilde(G) u^0||
$

$
(1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2 = (1 - delta) ||u'||^2 + delta ||u^r||^2 \
<= (1 - delta) (rho^(r/2) sqrt((1 - rho) / (2 - rho)) ||u^0||)^2 + delta (rho^(r/2) ||u^0||)^2 \
= (1 - delta) rho^r (1 - rho) / (2 - rho) lr(||u^0||)^2 + delta rho^r lr(||u^0||)^2 \
= rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta) lr(||u^0||)^2
$

Welche Werte kann $rho$ annehmen?


== Welche Werte kann $rho$ annehmen?

$
rho := (2 lambda^2) / (lambda^2 + 1) wide wide lambda := (|v^r|) / (||v^r||)
$
$
|v^r| <= ||v^r|| wide &==> wide 0 <= lambda <= 1
$

#place(center, dy: 20pt, {
	let xs = lq.linspace(0, 1, num: 200)
	lq.diagram(
		width: 325pt,
		height: 161pt,
		xaxis: (
			label: $lambda$,
		),
		yaxis: (
			label: $rho$,
		),
		lq.plot(
			xs,
			x => 2 * x * x / (x * x + 1),
			stroke: black + 4pt,
			mark: none,
		)
	)
})

==
Definiere
$
delta_q := max_(0 <= rho <= 1) rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)
$

$
(1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2 <= delta_q lr(||u^0||)^2
$
analog:
$
(1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2 <= delta_q ||hat(u)||^2
$

==

$
(hat(u), u^(2 r + 2)) <= sqrt((1 - delta) ||Q tilde(G) hat(u)||^2 + delta ||tilde(G) hat(u)||^2) dot sqrt((1 - delta) ||Q tilde(G) u^0||^2 + delta ||tilde(G) u^0||^2) \
<= sqrt(delta_q ||hat(u)||^2) dot sqrt(delta_q ||u^0||^2) = delta_q ||hat(u)|| ||u^0||
$


== Dualitätsargument
$
(hat(u), u^(2 r + 2)) <= delta_q ||hat(u)|| ||u^0|| wide forall hat(u) in S_h \
((hat(u), u^(2 r + 2))) / (||hat(u)||) = (hat(u) / (||hat(u)||), u^(2 r + 2)) <= delta_q ||u^0|| \
(hat(u) / (||hat(u)||), u^(2 r + 2)) <= underbrace(lr(||hat(u) / (||hat(u)||)||), = 1) dot ||u^(2 r + 2)|| \
||u^(2 r + 2)|| <= delta_q ||u^0||
$

== Ergebnis

$
||u^(2 r + 2)|| <= delta_q ||u^0||
$
mit
$
delta_q := max_(0 <= rho <= 1) rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta)
$
#pause
$
delta := cases(0 &wide& "gröbstes Level", (delta_(q - 1))^mu &wide& "feinere Level")
$

== Two-Grid-Konvergenz

$
delta = 0
$

$
delta_q = max_(0 <= rho <= 1) rho^r ((1 - delta) (1 - rho) / (2 - rho) + delta) \
	=^(delta = 0) max_(0 <= rho <= 1) rho^r (1 - rho) / (2 - rho)
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

#align(center + horizon,
	lq.diagram(
		width: 100%,
		height: 100%,
		..range(5).map(r => lq.plot(
			lq.linspace(1e-8, 1),
			x => calc.pow(x, r) * (1 - x) / (2 - x),
			mark: none,
			label: $r = #r$,
			stroke: 2pt,
		)),
		lq.scatter(
			maxima.map(m => m.at(0)),
			maxima.map(m => m.at(1)),
			color: black,
			size: 9pt,
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

== Multigrid-Konvergenz

$
delta = (delta_(q - 1))^mu \
delta_q = max_(0 <= rho <= 1) rho^r ((1 - (delta_(q - 1))^mu) (1 - rho) / (2 - rho) + (delta_(q - 1))^mu)
$
Annahme: $delta_q = delta_(q - 1) = delta_"MG"$ \
Löse die Gleichung
$
delta_"MG" = max_(0 <= rho <= 1) rho^r ((1 - (delta_"MG")^mu) (1 - rho) / (2 - rho) + (delta_"MG")^mu)
$

== Konvergenzraten

#figure(
	table(
		columns: 5,
		inset: 20pt,
		stroke: none,
		table.cell(align: right, $r ->$), table.vline(), $0$, $1$, $2$, $3$,
		table.hline(),
		[*Two Grid*], $0.5$, $0.172$, $0.114$, $0.086$,
		[*MG $mu = 1$ (V-cycle)*], $$, $1/2$, $1/3$, $1/4$,
		[*MG $mu = 2$ (W-cycle)*], $$, $0.187$, $0.120$, $0.087$
	)
)
