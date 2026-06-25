#set page(width: auto, height: auto, margin: 1cm)
#set text(font: "New Computer Modern", size: 11pt)

#import "@preview/cetz:0.3.4": canvas, draw

// ---------------------------------------------------------------------
// Konfiguration (gleiches Gitter wie triangulation.typ)
// ---------------------------------------------------------------------

#let n = 7
#let spacing = 2
#let node-radius = 0.55

// ---------------------------------------------------------------------
// Wiederverwendbare Zeichnung: Gitter + Knoten, deren Farbe über
// color-fn(i, j) -> color festgelegt wird
// ---------------------------------------------------------------------

#let mesh-canvas(color-fn) = align(center)[
  #canvas(length: 1cm, {
    import draw: *

    // Gitterlinien der regelmäßigen Triangulierung (Diagonalen gespiegelt)
    for j in range(0, n) {
      for i in range(0, n) {
        let p = (i * spacing, j * spacing)
        if i < n - 1 {
          line(p, ((i + 1) * spacing, j * spacing), stroke: gray + 0.6pt)
        }
        if j < n - 1 {
          line(p, (i * spacing, (j + 1) * spacing), stroke: gray + 0.6pt)
        }
        if i < n - 1 and j < n - 1 {
          if calc.rem(i + j, 2) == 0 {
            line(p, ((i + 1) * spacing, (j + 1) * spacing), stroke: gray + 0.6pt)
          } else {
            line((p.at(0) + spacing, p.at(1)), (p.at(0), p.at(1) + spacing), stroke: gray + 0.6pt)
          }
        }
      }
    }

    // Knoten als einfarbige Kreise
    for j in range(0, n) {
      for i in range(0, n) {
        let c = (i * spacing, j * spacing)
        circle(c, radius: node-radius, fill: color-fn(i, j), stroke: black + 0.5pt)
      }
    }
  })
]

// ---------------------------------------------------------------------
// Farbmuster für die drei Versionen
// ---------------------------------------------------------------------

// Version 1: alle Knoten blau
#let all-blue(i, j) = blue

// Version 2: Schachbrettmuster
// Spalte 0: Zeilen 0,2,4 blau; Spalte 1: Zeilen 1,3 blau; Spalte 2: Zeilen 0,2,4 blau; ...
#let checkerboard(i, j) = if calc.rem(i + j, 2) == 0 { blue } else { white }

// Version 3: Grobgitter
// nur jede zweite Spalte (0,2,4,...) hat in jeder zweiten Zeile (0,2,4,...) einen blauen Knoten
#let coarse-grid(i, j) = if calc.rem(i, 2) == 0 and calc.rem(j, 2) == 0 { blue } else { white }

// ---------------------------------------------------------------------
// Ausgabe
// ---------------------------------------------------------------------

#grid(
  columns: 3,
  gutter: 1.5cm,
  align(center)[*$Omega_h$* \ #mesh-canvas(all-blue)],
  align(center)[*$Omega_H$* \ #mesh-canvas(checkerboard)],
  align(center)[*$Omega_(2h)$* \ #mesh-canvas(coarse-grid)],
)

// ---------------------------------------------------------------------
// Einfarbiges Gitter mit dick/dünn-Kanten (für Folien, z.B. folien_anton.typ)
// ---------------------------------------------------------------------
// Im Unterschied zu mesh-canvas() oben markiert real-axis/real-diag, welche
// Kanten zur jeweiligen Norm "wirklich" gehören (dick) und welche nicht
// (dünn) -- damit lässt sich z.B. zeigen, dass Diagonalen nie zur Norm
// beitragen, oder dass nur die groben Kanten zur Seminorm beitragen.

#import "@preview/cetz:0.5.2" as cetz5
#import "lib.typ": W, H

#let thick-edge = 1.2pt
#let thin-edge = 0.25pt

#let single-color-grid(color-fn, real-axis: c => true, real-diag: (x, y) => true) = {
  import cetz5.draw: *
  for y in range(H + 1) {
    line((0, y), (W, y), stroke: if real-axis(y) { thick-edge } else { thin-edge })
  }
  for x in range(W + 1) {
    line((x, 0), (x, H), stroke: if real-axis(x) { thick-edge } else { thin-edge })
  }
  for x in range(W) {
    for y in range(H) {
      let stroke = if real-diag(x, y) { thick-edge } else { thin-edge }
      if calc.even(x + y) {
        line((x, y), (x + 1, y + 1), stroke: stroke)
      } else {
        line((x + 1, y), (x, y + 1), stroke: stroke)
      }
    }
  }
  for x in range(W + 1) {
    for y in range(H + 1) {
      circle((x, y), radius: 4pt, fill: color-fn(x, y), stroke: black + 0.4pt)
    }
  }
}

// Diagonale von Omega_(2h): abwechselnd nach links/rechts geneigt, je nach
// Position im 2x2-Block.
#let coarse-diag-real(x, y) = {
  let tmp = int(x / 2) + int(y / 2) + 1
  (calc.even(tmp) and calc.even(x + y)) or (calc.odd(tmp) and calc.odd(x + y))
}

// ---------------------------------------------------------------------
// Gitter für die Seminorm |v|
// ---------------------------------------------------------------------
// |v| summiert nur über die langen (2h-)Kanten von Omega_H; Diagonalen
// zählen nie. seminorm-grid() zeichnet genau das: das Omega_H-Schachbrett,
// bei dem nur diese 2h-Kanten dick (= "zählen für |v|") sind und alle
// anderen Kanten (kurze Kanten, Diagonalen) dünn bleiben.
#let seminorm-grid() = single-color-grid(
  checkerboard,
  real-axis: c => calc.even(c),
  real-diag: (x, y) => false,
)

#figure(
  cetz5.canvas({
    cetz5.draw.scale(1.5)
    seminorm-grid()
  }),
  caption: [Gitter für die Seminorm $|v|$ -- nur die dicken Kanten tragen bei.]
)
