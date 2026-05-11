#import "lib.typ": *

#import "@preview/typslides:1.3.3": *
#import "@preview/cetz:0.5.0"

#import calc: even, odd

#show: typslides.with(
	ratio: "16-9",
	theme: rgb("#0d3642"),
	font: "Fira Sans",
	font-size: 20pt,
	link-style: "color",
	show-progress: true,
)

#front-slide(
	title: "Convergence Rate of MG",
	subtitle: [With Gauss-Seidel Relaxation for the Poisson Equation],
	authors: "Oliver Portee",
	info: datetime.today().display("[day].[month].[year]"),
)

#table-of-contents()

#title-slide[
Grid
]


#let is-valid(x, y) = x + y >= 4 and x - y + 6 >= 2 and -y + 6 - x + 10 >= 2 and y - x + 10 >= 2
#let on-edge(x, y) = x + y == 4 or x - y + 6 == 2 or -y + 6 - x + 10 == 2 or y - x + 10 == 2

#slide[
#figure(
	cetz.canvas({
		import cetz.draw: *

		let dot(x, y) = {
			let fill = if even(x) and even(y) {
				green
			} else if even(x + y) {
				red
			} else {
				black
			}
			circle((x, y), radius: 4pt, stroke: none, fill: fill)
		}

	
		scale(2)

		for x in range(11) {
			for y in range(6) {
				if not is-valid(x, y) or not is-valid(x, y + 1) {
					continue
				}
				line((x, y), (x, y + 1), stroke: if even(x) { 3pt } else { 1pt })
			}
		}
		for y in range(7) {
			for x in range(10) {
				if not is-valid(x, y) or not is-valid(x + 1, y) {
					continue
				}
				line((x, y), (x + 1, y), stroke: if even(y) { 3pt } else { 1pt })
			}
		}
		for x in range(10) {
			for y in range(6) {
				if even(x + y) and is-valid(x, y) and is-valid(x + 1, y + 1) {
					line((x, y), (x + 1, y + 1), stroke: 3pt)
				}
				if odd(x + y) and is-valid(x + 1, y) and is-valid(x, y + 1) {
					line((x + 1, y), (x, y + 1), stroke: 3pt)
				}
			}
		}
		
		for x in range(11) {
			for y in range(7) {
				if not is-valid(x, y) {
					continue
				}
				dot(x, y)
			}
		}
	})
)
]
#let radius = 6pt
#let reddot = circle(radius: radius, stroke: none, fill: red)
#let greendot = circle(radius: radius, stroke: none, fill: green)
#let blackdot = circle(radius: radius, stroke: none, fill: black)

#slide(title: [Function Spaces])[
	$
		S_h &:= {u in op("span") {reddot, space greendot, space blackdot space "on fine grid"} | u = 0 "on" partial Omega} \
		S_H &:= {v in op("span") {reddot, space greendot space "on coarse grid"} | v = 0 "on" partial Omega} \
		T_h &:= {w in op("span") {blackdot space "on fine grid"} | w = 0 "on" partial Omega}
	$
	such that
	$
		S_h = S_H plus.o T_h \
		S_H subset S_h \
		T_h subset S_h
	$
]

// #let basis-function-image(divisions: 5, pixel-fn) = {
// 	let width = 10
// 	let height = 6
// 	let divisions = 5
// 	let W = divisions * width
// 	let H = divisions * height
// 	let d = (255,) * W * H
// 	let idx(x, y) = W * y + x
// 	let v(x, y) = is-valid(x, 6 - y)
// 	let e(vertices) = {
// 		let m(v1, v2) = {
// 			let (x1, y1) = v1
// 			let (x2, y2) = v2
// 			return ((x1 + x2) / 2, (y1 + y2) / 2)
// 		}
// 		let (a, b, c) = (vertices.at(0), vertices.at(1), vertices.at(2))
// 		let _edge(v) = is-edge(v.at(0), v.at(1))
// 		return _edge(m(a, b)) or _edge(m(b, c)) or _edge(m(a, c))
// 	}
// 	for yy in range(height) {
// 		for xx in range(width) {
// 			if even(xx + yy) {
// 				// SW of black
// 				if v(xx + 2/3, yy + 1/3) {
// 					let vertices = ((xx, yy), (xx + 1, yy), (xx + 1, yy + 1))
// 					for y in range(divisions) {
// 						for x in range(y, divisions) {
// 							let i = idx(divisions * xx + x, divisions * yy + y)
// 							d.at(i) = pixel-fn(xx, yy, x, y, divisions, "sw")	
// 						}
// 					}
// 				}
// 				// NE of black
// 				if v(xx + 1/3, yy + 2/3) {
// 					for y in range(divisions) {
// 						for x in range(y) {
// 							let i = idx(divisions * xx + x, divisions * yy + y)
// 							d.at(i) = pixel-fn(xx, yy, x, y, divisions, "ne")
// 						}
// 					}
// 				}
// 			} else {
// 				// SE of black
// 				if v(xx + 1/3, yy + 1/3) {
// 					for y in range(divisions) {
// 						for x in range(divisions - y) {
// 							let i = idx(divisions * xx + x, divisions * yy + y)
// 							d.at(i) = pixel-fn(xx, yy, x, y, divisions, "se")
// 						}
// 					}
// 				}
// 				// NW of black
// 				if v(xx + 2/3, yy + 2/3) {
// 					for y in range(divisions) {
// 						for x in range(divisions - y, divisions) {
// 							let i = idx(divisions * xx + x, divisions * yy + y)
// 							d.at(i) = pixel-fn(xx, yy, x, y, divisions, "nw")
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}
// 	image(
// 		bytes(d),
// 		format: (
// 			encoding: "luma8",
// 			width: divisions * width,
// 			height: divisions * height,
// 		),
// 		width: 37%,
// 	)
// }


#let g(v) = calc.pow(v, 0.8)

// #let pixel_Th(xx, yy, x, y, divisions, dir) = {
// 	if dir == "sw" {
// 		(x, y) = (divisions - x, y)
// 	} else if dir == "se" {
// 		(x, y) = (x, y)
// 	} else if dir == "nw" {
// 		(x, y) = (divisions - x, divisions - y)
// 	} else if dir == "ne" {
// 		(x, y) = (x, divisions - y)
// 	}
// 	int(g((x + y) / divisions) * 255)
// }

// #let pixel_SH-red(xx, yy, x, y, divisions, dir) = {
// 	let a = 0
// 	if (dir == "sw" or dir == "se") and even(yy) {
// 		a = divisions - y
// 	} else if (dir == "nw" or dir == "ne") and odd(yy) {
// 		a = y
// 	} else if (dir == "nw" or dir == "sw") and odd(xx) {
// 		a = x
// 	} else if (dir == "ne" or dir == "se") and even(xx) {
// 		a = divisions - x
// 	}
// 	int(g(a / divisions) * 255)	
// }


// #basis-function-image(pixel_Th)

// #basis-function-image(pixel_SH-red)

#let width = 10
#let height = 6
#let divisions = 20
#let W = divisions * width
#let H = divisions * height
#let idx(x, y) = W * y + x
#let v(x, y) = is-valid(x, 6 - y)
#let image-width = 80%

#let T_h = {
	let black-hat(d, xx, yy) = {
		if not v(xx, yy) {
			return d
		}
		for y in range(divisions) {
			for x in range(divisions - y) {
				let (a, b) = (x, y)
				let i = idx(divisions * xx + a, divisions * yy + b)
				d.at(i) = int(g((x + y) / divisions) * 255)
				let (a, b) = (divisions - x, y)
				let i = idx(divisions * (xx - 1) + a, divisions * yy + b)
				d.at(i) = int(g((x + y) / divisions) * 255)
				let (a ,b) = (x, divisions - y)
				let i = idx(divisions * xx + a, divisions * (yy - 1) + b)
				d.at(i) = int(g((x + y) / divisions) * 255)
				let (a ,b) = (divisions - x, divisions - y)
				let i = idx(divisions * (xx - 1) + a, divisions * (yy - 1) + b)
				d.at(i) = int(g((x + y) / divisions) * 255)
				
			}
		}
		return d
	}

	let d = (255,) * W * H
	for y in range(1, height) {
		for x in range(calc.rem(y + 1, 2), width, step: 2) {
			d = black-hat(d, x, y)		
		}
	}

	image(
		bytes(d),
		format: (
			encoding: "luma8",
			width: divisions * width,
			height: divisions * height,
		),
		width: image-width,
	)
}



#let S_H-red = {
	let red-hat(d, xx, yy) = {
		if not v(xx, yy) or on-edge(xx, height - yy) {
			return d
		}
		for y in range(-divisions, divisions) {
			for x in range(calc.abs(y), divisions) {
				let (a, b) = (x, y)
				let i = idx(divisions * xx + a, divisions * yy + b)
				d.at(i) = int(g(x / divisions) * 255)
				let (a, b) = (divisions - x, y)
				let i = idx(divisions * (xx - 1) + a, divisions * yy + b)
				d.at(i) = int(g(x / divisions) * 255)
				let (a, b) = (y, x)
				let i = idx(divisions * xx + a, divisions * yy + b)
				d.at(i) = int(g(x / divisions) * 255)
				let (a, b) = (y, divisions - x)
				let i = idx(divisions * xx + a, divisions * (yy - 1) + b)
				d.at(i) = int(g(x / divisions) * 255)
				
			}
		}
		return d
	}

	let d = (255,) * W * H
	for y in range(1, height, step: 2) {
		for x in range(1, width, step: 2) {
			d = red-hat(d, x, y)		
		}
	}

	image(
		bytes(d),
		format: (
			encoding: "luma8",
			width: divisions * width,
			height: divisions * height,
		),
		width: image-width,
	)
}



#let S_H-green-1 = {
	let D = 2 * divisions

	let green-hat(d, xx, yy) = {
		if not v(xx, yy) or on-edge(xx, height - yy) {
			return d
		}
		for y in range(D) {
			for x in range(D - y) {
				let (a, b) = (x, y)
				let i = idx(divisions * xx + a, divisions * yy + b)
				d.at(i) = int(g((x + y) / D) * 255)
				let (a, b) = (D - x, y)
				let i = idx(divisions * (xx - 2) + a, divisions * yy + b)
				d.at(i) = int(g((x + y) / D) * 255)
				let (a, b) = (x, D - y)
				let i = idx(divisions * xx + a, divisions * (yy - 2) + b)
				d.at(i) = int(g((x + y) / D) * 255)
				let (a, b) = (D - x, D - y)
				let i = idx(divisions * (xx - 2) + a, divisions * (yy - 2) + b)
				d.at(i) = int(g((x + y) / D) * 255)
				
			}
		}
		return d
	}

	let d = (255,) * W * H
	for y in range(2, height, step: 2) {
		let xstart = calc.rem(int(y / 2), 2) * 2
		for x in range(xstart, width, step: 4) {
			d = green-hat(d, x, y)		
		}
	}

	image(
		bytes(d),
		format: (
			encoding: "luma8",
			width: divisions * width,
			height: divisions * height,
		),
		width: image-width,
	)
}

#let S_H-green-2 = {
	let D = 2 * divisions

	let green-hat(d, xx, yy) = {
		if not v(xx, yy) or on-edge(xx, height - yy) {
			return d
		}
		for y in range(D) {
			for x in range(D - y) {
				let (a, b) = (x, y)
				let i = idx(divisions * xx + a, divisions * yy + b)
				d.at(i) = int(g((x + y) / D) * 255)
				let (a, b) = (D - x, y)
				let i = idx(divisions * (xx - 2) + a, divisions * yy + b)
				d.at(i) = int(g((x + y) / D) * 255)
				let (a, b) = (x, D - y)
				let i = idx(divisions * xx + a, divisions * (yy - 2) + b)
				d.at(i) = int(g((x + y) / D) * 255)
				let (a, b) = (D - x, D - y)
				let i = idx(divisions * (xx - 2) + a, divisions * (yy - 2) + b)
				d.at(i) = int(g((x + y) / D) * 255)
				
			}
		}
		return d
	}

	let d = (255,) * W * H
	for y in range(2, height, step: 2) {
		let xstart = calc.rem(int(y / 2) + 1, 2) * 2
		for x in range(xstart, width, step: 4) {
			d = green-hat(d, x, y)		
		}
	}

	image(
		bytes(d),
		format: (
			encoding: "luma8",
			width: divisions * width,
			height: divisions * height,
		),
		width: image-width,
	)
}

#grid(
	columns: (1fr, 1fr),
	align: center + horizon,
	row-gutter: 10pt,
	T_h, S_H-red,
	S_H-green-1, S_H-green-2,
)

system of equations
$
	u_i = 1/4 attach(sum, b: j, br: h, tr: ') u_j + b_i wide p_i in Omega_h
$

#slide(title: [Splitting the GS Step])[
	$
		(G_h^I u)_i &= cases(u_i\, &wide& p_i "is" greendot "or" reddot, 1/4 sum_("neighbor" j) u_j + b_i\, &wide& "otherwise") \
		(G_h^(I I) u)_i &= cases(1/4 sum_("neighbor" j) u_j + b_i\, &wide& p_i "is" greendot "or" reddot, u_i &wide& "otherwise")
	$
]

#slide[
#align(center,
	cetz.canvas({
		import cetz.draw: *

		scale(y: -100%)

		h-grid()
		content((W / 2, H + 1), $Omega_h$)
		translate(x: 11)
		H-grid()
		content((W / 2, H + 1), $Omega_H$)
		translate(x: -11, y: 7)
		h2-grid()
		content((W / 2, H + 1), $Omega_(2 h)$)
		translate(x: 11)
		content((W / 2, H / 2), $H = sqrt(2) h$)
	})
)
]

#slide[
#align(center,
	cetz.canvas({
		import cetz.draw: *

		scale(y: -100%)

		h-grid()
		content((W / 2, H + 1), $u in S_h$)
		content((8, H / 2), $=$)
		translate(x: 10)
		H-grid()
		content((W / 2, H + 1), $v in S_H$)
		content((8, H / 2), $plus.o$)
		translate(x: 10)
		h-grid(hide-green-red: true)
		content((W / 2, H + 1), $w in T_h$)
	})
)
]

#slide[
$
integral_I v_xi w_xi + v_eta w_eta = thin ?
$
$
integral_I v_xi w_xi = integral_(I_1) v_xi w_xi + integral_(I_2) v_xi w_xi = v_xi integral_(I_1) w_xi + v_xi integral_(I_2) w_xi
$
since $lr(w_xi|)_(I_1) = - lr(w_xi|)_(I_2)$ and $lr(v_xi|)_I = "const"$
$
integral_I v_xi w_xi = 0
$
]

#slide[

]
