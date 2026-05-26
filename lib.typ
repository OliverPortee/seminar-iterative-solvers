#import "@preview/cetz:0.5.0"
#import "@preview/physica:0.9.8": *

#import calc: even, odd

#let W = 6
#let H = 4

#let grid-point(x, y) = {
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
					grid-point(x, y)
				}
			} else {
				grid-point(x, y)
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
				grid-point(x, y)
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
			grid-point(x, y)
		}
	}
}

