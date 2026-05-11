#import "lib.typ": *

#set heading(numbering: "1.")

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