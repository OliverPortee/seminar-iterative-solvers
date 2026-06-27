#import "@preview/touying:0.7.4": *
#import "@preview/cetz:0.5.2"
#import "@preview/lilaq:0.6.0" as lq

#import calc: even, odd
#import themes.simple: *

#import "lib.typ": W, H, draw-grid, draw-coarse-grid, draw-2h-grid, draw-th-grid, fine-color, coarse-color

#show: simple-theme

#let hi(content) = box(width: 100%, inset: 0.8em, stroke: orange, content)

#let Gbox(content, color) = box(fill: color.transparentize(70%), inset: (top: 7pt, bottom: 4pt, left: 3pt, right: 2pt), radius: 5pt, stroke: color, baseline: 15%, content)
#let G1(content) = Gbox(content, fine-color)
#let G2(content) = Gbox(content, coarse-color)

// einfarbige Gitter wie im Handout (statt der bunten lib.typ-Gitter)
// single-color-grid/seminorm-grid liegen zentral in triangulation-single-color.typ
#import "triangulation-single-color.typ": all-blue, single-color-grid, seminorm-grid

#title-slide[
#v(1fr)
= MG Konvergenzrate
Teil 1: Grundlagen

#v(1fr)
Anton Wiede
]

= Problemstellung

== Die Poisson-Gleichung

$
-laplace u = f &wide& "in" Omega \
u = 0 &wide& "auf" partial Omega
$

#pause



*FEM*: bringe Poisson Gleichung zuerst in die schwache Form

$
	integral_Omega nabla u dot nabla v = integral_Omega f dot v wide forall v in S_h
$

== Diskretisierung

- $Omega$ wird trianguliert
- $u$ wird stückweise linear auf den Gitterpunkten approximiert
- statt der PDE lösen wir ein (großes, dünnes) lineares Gleichungssystem

#hi[Idee von Multigrid: löse das Problem nicht nur auf einem Gitter, sondern auf einer ganzen Hierarchie von immer gröberen Gittern.]

= Gitterhierarchie

== Drei Gitter

#figure(
	cetz.canvas({
		import cetz.draw: content, translate, scale
		scale(1)
		draw-grid()
		translate(x: W + 2)
		draw-coarse-grid()
		translate(x: W + 2)
		draw-2h-grid()

	})
)


- #G1[$Omega_h$]: Kantenlängen $h$ und $sqrt(2) h$ — das feine Gitter
- #G2[$Omega_H$]: Kantenlängen $sqrt(2) h$ und $2h$ — jeder zweite Punkt
- $Omega_(2h)$: Kantenlängen $2h$ und $2 sqrt(2) h$ — noch gröber


== Warum mehrere Gitter?

#hi[
Fehler mit kurzer Wellenlänge (hochfrequent) werden durch lokale Relaxation (Gauß-Seidel) schnell geglättet.

Fehler mit langer Wellenlänge (niederfrequent) sieht man auf dem feinen Gitter kaum — aber auf einem groben Gitter sind sie wieder hochfrequent!
]

#G1[fein] glättet die kurzen Wellen #h(1fr) #G2[grob] übernimmt die langen Wellen

= Raumzerlegung

== Finite-Elemente-Räume

$
u in S_h &: "Finite-Elemente-Raum auf" Omega_h \
v in S_H &: "Finite-Elemente-Raum auf" Omega_H \
w in T_h &:= {w in S_h | w(p) = 0 wide forall p in Omega_H}
$

*Intuition:* $T_h$ enthält genau die "Details", die nur auf dem feinen Gitter sichtbar sind — auf jedem groben Gitterpunkt verschwinden sie.

$
S_h &= S_H plus.o T_h wide wide wide  u = v+w
$

== Finite-Elemente-Raum

#figure(
	cetz.canvas({
		import cetz.draw: content, translate, scale
		scale(1.2)
		draw-grid()
		content((7, 2), $=$)
		content((W / 2, -0.5), $u in S_h$)
		translate(x: 8)
		draw-coarse-grid()
		content((7, 2), $plus.o$)
		content((W / 2, -0.5), $v in S_H$)
		translate(x: 8)
		draw-th-grid()
		content((W / 2, -0.5), $w in T_h$)
	}),
	caption: [$S_h$ (fein) zerfällt in $S_H$ (grob, gelb) und $T_h$ (blau, nur dazwischen).]
)


= Energienorm

== Die Bilinearform $a(u, v)$

$
a(u, v) := integral_Omega nabla u dot nabla v = integral_Omega (u_xi v_xi + u_eta v_eta) thin dif xi thin dif eta
$

*Intuition:* $a(dot, dot)$ ist ein Skalarprodukt für Funktionen — so wie $x dot y$ im $RR^n$, nur dass statt "Richtung" die *Steigung* (Gradient) verglichen wird.

#hi[a(u,v) entspricht der linken Seite in der schwachen Form]

#pause


== Auswertung der Energienorm

$
||u|| := sqrt(a(u, u)) wide "(Energienorm, die zugehörige" "Länge")
$
#figure(
	cetz.canvas({
		import cetz.draw: *
		scale(1.2)
		single-color-grid(all-blue, real-diag: (x, y) => false)
	})
)

$
||u||^2 = a(u,u) = sum_vec(i\, j in Omega_h, d(i,j) = h, delim: #none) (u_i - u_j)^2
$


== Die Seminorm

#grid(
	columns: (auto, 1fr),
	align: horizon,
	column-gutter: 1.5em,
	figure(
		cetz.canvas({
			import cetz.draw: *
			scale(1.1)
			seminorm-grid()
		})
	),
	$
	|u|^2 := sum_vec(i\, j in Omega_H, d(i,j) = 2h, delim: #none) 1/2 (u_i - u_j)^2
	$
)

#hi[Das ist eine *Seminorm*: $|u| = 0$ ist möglich, obwohl $u != 0$ — nämlich wenn $u$ nur auf den schwarzen (feinen) Punkten variiert.]

== Ziel: Konvergenzrate bestimmen

Wir wollen eine Konstante finden sodass gilt:

$
  ||u^(k+1)|| <= C * ||u||
$

== Konvergenzrate Gauss-Seidel

Hierfür betrachten wir zuerst die Gauss-Seidel Relaxation auf dem feinsten Gitter. Hier wollen wir zeigen, dass

$
(1) wide |u^k| <= sqrt(rho) dot ||u^k|| wide "(Der grobe Fehler ist kleiner als der feine)" \

(2) wide ||u^(k+1)|| <= |u^k| wide "(der feine Fehler" \ " im nächsten Schritt ist kleiner als der Grobe Fehler im vorherigen)"
$

$
==> wide ||u^(k+1)|| <= sqrt(rho) dot ||u^k||
$

== Vergleich $|u| <= ||u||$

#pause

- $|u|$ zählt die Hälfte der Kanten von  $||u||$
#pause
- jede grobe Kante lässt sich durch zwei feine Kanten "ersetzen" 
#pause
$
1/2 (c-a)^2 <=(c-b)^2 + (b-a)^2 ==>  |u| <= ||u|| wide forall u in S_h
$

#pause

= Obere Schranke für $(|u|) / (||u||)$

== Herleitung der Schranke

#grid(
	columns: (50%, 1fr),
	align: horizon,
	column-gutter: 1.5em,
	figure(
		cetz.canvas({
			import cetz.draw: line, content, scale, circle

			scale(3)

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
			circle(p1, radius: r, stroke: none, fill: yellow)
			circle(p2, radius: r, stroke: none, fill: yellow)
			circle(p3, radius: r, stroke: none, fill: yellow)
			circle(p4, radius: r, stroke: none, fill: yellow)
			circle(p5, radius: r, stroke: none, fill: blue)

			content(p1, $1$, anchor: "east", padding: 0.25)
			content(p3, $3$, anchor: "east", padding: 0.25)
			content(p2, $2$, anchor: "north", padding: 0.25)
			content(p4, $4$, anchor: "north", padding: 0.25)
			content(p5, $5$, anchor: "north-east", padding: 0.2)

			content((0.5, 0), $h$, anchor: "south", padding: 0.1)
			content((0, 0.5), $h$, anchor: "west", padding: 0.1)
		})
	),
	[
		- $Omega$ wird in Dreieckspaare $"I" = triangle 1 2 4$, $"II" = triangle 4 2 3$ mit gemeinsamer Kante $overline(4 2)$ zerlegt
	],
)

== Herleitung der Schranke
- $integral_"I" v_xi w_xi = 0$, da $w_xi$ auf den beiden Hälften von $"I"$ entgegengesetztes Vorzeichen
- $v_xi$ ist konstant
- quadratische Ergänzung von $integral_"I" v_eta w_eta$ liefert folgendes Ergebnis

$
a(v,w) <= 1/2 (||v||^2 - |v|^2) + 1/4 ||w||^2
$

== Herleitung der Schranke

da $a(v,w)$ linear in $w$ und die rechte Seite quadratisch in $w$ ist:
$
a(v,w) <= 1/2 (||v||^2 - |v|^2) + 1/4 ||w||^2
$
$
==> |a(v, w)| <= sqrt(1/2 thin (||v||^2 - |v|^2)) ||w|| = sqrt(1/2 thin (1 - underbrace((|v|^2) / (||v||^2), lambda^2))) ||v||||w||.
$

== Geometrische Interpretation

#figure(
	cetz.canvas({
		import cetz.draw: line, content, arc, scale
		scale(1.3)
		let a = (0, 0)
		let b = (6, 0)
		let c = (4, 2)
		line(a, b, mark: (end: ">"), name: "v", stroke: 2pt + coarse-color)
		line(b, c, mark: (end: ">"), name: "w", stroke: 2pt + fine-color)
		line(a, c, mark: (end: ">"), name: "u", stroke: 2pt)
		content("v", $v$, padding: 0.1, anchor: "north")
		content("u", $u$, padding: 0.1, anchor: "south-east")
		content("w", $w$, padding: 0.1, anchor: "south-west")
		arc(b, anchor: "origin", radius: 1, start: 135deg, stop: 180deg)
		content((5.4, 0.28), $alpha$)
	})
)

$
cos(alpha) = a(v,w) / (||v|| ||w||)
$

Mit dem Skalarprodukt $a(dot,dot)$ verhält sich $S_h$ wie ein gewöhnlicher euklidischer Raum.

== Das Lemma

$
|a(v,w)| <= sqrt(1/2 (1 - lambda^2)) thin ||v|| thin ||w|| wide forall v in S_H, w in T_h
$

äquivalent (mit der Winkel-Interpretation):
$
cos^2(alpha) <= 1/2 (1 - lambda^2) \
==> sin^2 alpha >= 1/2 (1 + lambda^2)
$

== Ergebnis

Generell gilt $||u|| >= ||v|| sin(alpha)$

$
==> |u| = |v| = lambda||v|| <= lambda / sin(alpha)||u|| <= lambda / (sqrt(1/2 (1 + lambda^2)))||u||
$
$
==> (|u|) / (||u||) <= lambda / sqrt(1/2 (1 + lambda^2)) = sqrt(lambda^2 / (1/2 (1 + lambda^2))) = sqrt((2 lambda^2) / (1 + lambda^2)) =: sqrt(rho)
$


== $rho$ als Konvergenzparameter

$
rho := (2 lambda^2) / (1 + lambda^2)
$

#figure({
	let xs = lq.linspace(0, 1, num: 200)
	lq.diagram(
		width: 320pt,
		height: 160pt,
		xaxis: (label: $lambda$),
		yaxis: (label: $rho$),
		lq.plot(xs, x => 2 * x * x / (x * x + 1), stroke: black + 3pt, mark: none),
	)
})

== Zwischenstand: die Konvergenzrate

$
&(1) wide |u^k| <= sqrt(rho) dot ||u^k|| wide & "(gerade gezeigt)" \
&(2) wide ||u^(k+1)|| <= |u^k| wide & "(folgt aus Gauß-Seidel)"
$

$
==> wide ||u^(k+1)|| <= sqrt(rho) dot ||u^k||
$



= Gauß-Seidel-Relaxation

== Der Stencil

Für $f = 0$ ergibt die finite-Elemente-Diskretisierung an jedem Punkt
$
u_i = 1/4 sum_(j in N(i)) u_j,
wide
mat(0, -1, 0; -1, 4, -1; 0, -1, 0)
$

*Intuition:* die Lösung der homogenen Gleichung ist an jedem Punkt der Mittelwert der vier Nachbarn — genau das berechnet ein Gauß-Seidel-Schritt.

== Die zwei Halbschritte

$
(G_h^"I" u)_i &= cases(
	u_i &wide p_i in Omega_H,
	1/4 sum_(j in N(i)) u_j &wide p_i in Omega_h \\ Omega_H,
) \
(G_h^"II" u)_i &= cases(
	1/4 sum_(j in N(i)) u_j &wide p_i in Omega_H,
	u_i &wide p_i in Omega_h \\ Omega_H,
)
$
$G_h^"II"$ glättet nur auf den groben Gitterpunkten \
$G_h^"I"$ glättet auf allen anderen

$N(i)$ ist die Menge aller benachbarten Punkte mit Abstand $h$ von $p_i$.

== Konvergenz

*Lemma * $Omega$ konvex. Falls $u = G_h^"I" u$, dann
$
 ||G_h^"II" u|| <= |u|.
$

#hi[Intuition: Wenn die feinen Punkte bereits geglättet sind, reicht eine Glättung der groben Punkte um die Energienorm auf die Seminorm zu senken.]



#let diamond-fig() = cetz.canvas({
	import cetz.draw: *
	scale(1.5)
	let p0 = (0, 0)
	let p1 = (0, -1)
	let p2 = (1, 0)
	let p3 = (0, 1)
	let p4 = (-1, 0)
	line(p1, p2, p3, p4, close: true)
	line(p0, p1, stroke: (dash: "dashed"))
	line(p0, p2, stroke: (dash: "dashed"))
	line(p0, p3, stroke: (dash: "dashed"))
	line(p0, p4, stroke: (dash: "dashed"))
	circle(p0, radius: 5pt, fill: coarse-color, stroke: black)
	circle(p1, radius: 5pt, fill: fine-color, stroke: black)
	circle(p2, radius: 5pt, fill: fine-color, stroke: black)
	circle(p3, radius: 5pt, fill: fine-color, stroke: black)
	circle(p4, radius: 5pt, fill: fine-color, stroke: black)
	content(p0, $0$, anchor: "south", padding: 0.2)
	content(p1, $1$, anchor: "east", padding: 0.3)
	content(p2, $2$, anchor: "west", padding: 0.3)
	content(p3, $3$, anchor: "west", padding: 0.3)
	content(p4, $4$, anchor: "east", padding: 0.3)
})

== Schritt 1: Lokale Identität

#grid(
	columns: (50%, 1fr),
	align: horizon,
	column-gutter: 1.5em,
	figure(diamond-fig()),
	[
		- #G2[$0$]: Punkt in $Omega_H$, gestrichelt: Abstand $h$
		- #G1[$1,2,3,4$]: seine vier Nachbarn in $Omega_h without Omega_H$
	],
)



Sei $macron(u) := G_h^"II" u$ . Da $G_h^"II"$ den Punkt $0$ auf den Mittelwert seiner Nachbarn setzt,
$
macron(u)_0 = 1/4 sum_(i=1)^4 u_i = 1/4 sum_(i=1)^4 macron(u)_i
$


== Schritt 2: Zerlegung der Energienorm $Omega_H$


Wie bei der Energienorm zerlegt nach Gitterebene:
$
||macron(u)||^2 = sum_(i in Omega_H) sum_(j in N(i)) (macron(u)_i - macron(u)_j)^2
$

== Schritt 2: Zerlegung der Energienorm $Omega_H$

Umformen liefert:
$
||G_h^"II" u||^2 <= sum_(k,l in Omega_h without Omega_H, \ d(k,l)=H) (u_k - u_l)^2 overbrace(<=, u = G_h^"I"u) 1/4 sum_(m,n in Omega_H, \ d(k,l)=2H) (u_m - u_n)^2 \ = 1/2 sum_(n,m in Omega_H, \ d(n,m)=2h) (u_m - u_n)^2 = |u|
$


== Finaler Schritt

Sei $Omega$ konvex. Dann gilt für *jedes* $u in S_h$
$
||G_h^"II" dot G_h^"I" u|| <= |u|.
$


*Beweis:* Setze $u' := G_h^"I" u$. Dann ist $u' = G_h^"I" u'$ (Idempotenz), also greift das Lemma auf $u'$:
$
||G_h^"II" u'|| <= |u'|
$


Außerdem ändert $G_h^"I"$ nur den $w$-Anteil ($T_h$-Teil), nicht den $v$-Anteil — und $|dot|$ hängt nur vom $v$-Anteil ab (vgl. $|u| = |v|$ aus der Zerlegung $u = v + w$):
$
|u'| = |G_h^"I" u| = |u|
$

$
==> ||G_h^"II" G_h^"I" u|| <= |u| wide forall u in S_h. wide qed
$