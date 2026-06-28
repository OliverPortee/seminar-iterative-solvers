#import "@preview/cetz:0.5.2"
#import "@preview/physica:0.9.8": *
#import "@preview/lilaq:0.6.0" as lq
#import "@preview/tiptoe:0.4.0"
#import "@preview/mannot:0.3.3": markhl

#import calc: even, odd

#let W = 6
#let H = 4

// Klassisches dreifarbiges Schema (oliv/rot/schwarz), gebraucht von
// h-grid/H-grid/h2-grid unten (handout.typ).
#let grid-point-classic(x, y) = {
	let fill = if even(x) and even(y) {
		olive
	} else if odd(x) and odd(y) {
		red
	} else {
		black
	}
	cetz.draw.circle((x, y), fill: fill, stroke: none, radius: 4pt)
}

#let h-grid(hide-green-red: false) = {
	for y in range(H + 1) {
		cetz.draw.line((0, y), (W, y))
	}
	for x in range(W + 1) {
		cetz.draw.line((x, 0), (x, H))
	}
	for x in range(W) {
		for y in range(H) {
			if even(x + y) {
				cetz.draw.line((x, y), (x + 1, y + 1))
			} else {
				cetz.draw.line((x + 1, y), (x, y + 1))
			}
		}
	}
	for x in range(W + 1) {
		for y in range(H + 1) {
			if even(x + y) {
				if not hide-green-red {
					grid-point-classic(x, y)
				}
			} else {
				grid-point-classic(x, y)
			}
		}
	}
}

#let H-grid() = {
	for y in range(0, H + 1, step: 2) {
		cetz.draw.line((0, y), (W, y))
	}
	for x in range(0, W + 1, step: 2) {
		cetz.draw.line((x, 0), (x, H))
	}
	for x in range(W) {
		for y in range(H) {
			if even(x + y) {
				cetz.draw.line((x, y), (x + 1, y + 1))
			} else {
				cetz.draw.line((x + 1, y), (x, y + 1))
			}
		}
	}
	for x in range(W + 1) {
		for y in range(H + 1) {
			if even(x + y) {
				grid-point-classic(x, y)
			}
		}
	}
}

#let h2-grid() = {
	for y in range(0, H + 1, step: 2) {
		cetz.draw.line((0, y), (W, y))
	}
	for x in range(0, W + 1, step: 2) {
		cetz.draw.line((x, 0), (x, H))
	}
	for x in range(W) {
		for y in range(H) {
			let tmp = int(x / 2) + int(y / 2) + 1 // remove "+ 1" to use the other connections
			if even(tmp) and even(x + y) {
				cetz.draw.line((x, y), (x + 1, y + 1))
			} else if odd(tmp) and odd(x + y) {
				cetz.draw.line((x + 1, y), (x, y + 1))
			}
		}
	}
	for x in range(0, W + 1, step: 2) {
		for y in range(0, H + 1, step: 2) {
			grid-point-classic(x, y)
		}
	}
}

// Neueres zweifarbiges Schema (blau/gelb = fein/grob), gebraucht von
// draw-grid/draw-coarse-grid/draw-2h-grid/draw-th-grid unten
// (handout_anton.typ, folien_anton.typ).
#let fine-color = lq.color.map.petroff10.at(0)
#let coarse-color = lq.color.map.petroff10.at(1)

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

