using SymPy

#==
bookerquartic
==#
@vars θ M11 M12 M13 M21 M22 M23 M31 M32 M33 dC² dS
S = sin(θ)
C = cos(θ)
C² = C^2
B4 = 1 + M33
B3 = S*(M13 + M31)
B2 = -(C² + M33)*(1 + M11) + M13*M31 -
     (1 + M33)*(C² + M22) + M23*M32
B1 = S*(M12*M23 + M21*M32 -
     (C² + M22)*(M13 + M31))
B0 = (1 + M11)*(C² + M22)*(C² + M33) +
     M12*M23*M31 + M13*M21*M32 -
     M13*(C² + M22)*M31 -
     (1 + M11)*M23*M32 -
     M12*M21*(C² + M33)

# dS = cos(θ)
# dC² = 2*sin(θ)*cos(θ)

diff(C², θ)
dB4 = diff(B4, θ)
dB3 = diff(B3, θ)
dB2 = diff(B2, θ).subs(2*sin(θ)*cos(θ), -dC²)
dB1 = diff(B1, θ).subs(2*sin(θ)*cos(θ), -dC²).subs(cos(θ), dS)
dB0 = simplify(diff(B0, θ).subs(2*sin(θ)*cos(θ), -dC²))

@vars q
tmp = diff(B4*q^4 + B3*q^3 + B2*q^2 + B1*q + B0, θ).subs(2*sin(θ)*cos(θ), -dC²)
tmp = tmp.subs((M13 + M31)*cos(θ), Sym("dB3")).subs(-dC²*(M11+1)+dC²*(-M33-1), Sym("dB2"))
tmp = tmp.subs(dC²*(-M13 - M31)*sin(θ) + cos(θ)*(M12*M23 + M21*M32 - (M13 + M31)*(M12 + C²)), Sym("dB1"))
tmp.subs(dC²*(-M12*M21 - M13*M31 + (M11 + 1)*(M22 + C²) +
    (M11 + 1) * (M33 + C²)), Sym("dB0"))


@vars θ
q1 = sympy.Function("q1")(θ)
q2 = sympy.Function("q2")(θ)
den = sympy.Function("den")(θ)
C = cos(θ)

T1 = sympy.Function("T1")(θ)
T2 = sympy.Function("T2")(θ)
P1 = sympy.Function("P1")(θ)
P2 = sympy.Function("P2")(θ)

diff(-2*C*(q1-q2)*den, θ)  # R21

diff(-2*C*(T1*P2 - T2*P1)*den, θ)  # R12

# Ground R

@vars θ ng²

C = cos(θ)
S = sin(θ)

tmp1 = C*ng²
tmp2 = sqrt(ng² - S^2)

# diff((C*ng^2 - sqrt(ng^2 - S^2))/(C*ng^2 + sqrt(ng^2 - S^2)), θ)

diff((tmp1 - tmp2)/(tmp1 + tmp2), θ)  # rg11

diff((C - tmp2)/(C + tmp2), θ)  # rg22  factor() for simplification


# height gain function
@vars z d a h10 h20 k

h1 = sympy.Function("h1")(z)
h2 = sympy.Function("h2")(z)
F1 = sympy.Function("F1")(z)
F2 = sympy.Function("F2")(z)

f1 = exp(z/a)*(F1*h1 + F2*h2)/(F1*h10 + F2*h20)
1/(im*k)*diff(f1, z)

# If dF/dz = 0
@vars z d a h10 h20 k F1 F2

h1 = sympy.Function("h1")(z)
h2 = sympy.Function("h2")(z)

f1 = exp(z/a)*(F1*h1 + F2*h2)/(F1*h10 + F2*h20)
1/(im*k)*diff(f1, z)


# = -im/(a*k)*exp(z/a)*(F1*h1(z) + F2*h2(z) + a*(F1*h1p(z) + F2*h2p(z)))/(F1*h10 + F2h20)

f₁ = exp(z/Rₑ)*(F₁*h₁z + F₂*h₂z)/(F₁*h₁0 + F₂*h₂0)


@vars k z Re F₁ F₂

h1 = sympy.Function("h1")(z)
h2 = sympy.Function("h2")(z)

f₁ = exp(z/Re)*(F₁*h1 + F₂*h2)

1/(im*k)*diff(f₁, z)


#
# dRdz
#
@vars k θ

W1_11 = sympy.Function("W1_11")(θ)
W1_12 = sympy.Function("W1_12")(θ)
W1_21 = sympy.Function("W1_21")(θ)
W1_22 = sympy.Function("W1_22")(θ)

W2_11 = sympy.Function("W2_11")(θ)
W2_12 = sympy.Function("W2_12")(θ)
W2_21 = sympy.Function("W2_21")(θ)
W2_22 = sympy.Function("W2_22")(θ)

W3_11 = sympy.Function("W3_11")(θ)
W3_12 = sympy.Function("W3_12")(θ)
W3_21 = sympy.Function("W3_21")(θ)
W3_22 = sympy.Function("W3_22")(θ)

W4_11 = sympy.Function("W4_11")(θ)
W4_12 = sympy.Function("W4_12")(θ)
W4_21 = sympy.Function("W4_21")(θ)
W4_22 = sympy.Function("W4_22")(θ)

R11 = sympy.Function("R11")(θ)
R21 = sympy.Function("R21")(θ)
R12 = sympy.Function("R12")(θ)
R22 = sympy.Function("R22")(θ)

R = [R11 R12;
     R21 R22]

W1 = [W1_11 W1_12;
      W1_21 W1_22]
W2 = [W2_11 W2_12;
      W2_21 W2_22]
W3 = [W3_11 W3_12;
      W3_21 W3_22]
W4 = [W4_11 W4_12;
      W4_21 W4_22]

@vars k θ

R = sympy.Function("R")(θ)
W1 = sympy.Function("W1")(θ)
W2 = sympy.Function("W2")(θ)
W3 = sympy.Function("W3")(θ)
W4 = sympy.Function("W4")(θ)

dz = -im/2*k*(W2 + W4*R - R*W1 - R*W3*R)

diff(dz, θ)


dθdz = i/2*k*(-R*dW3*R + 2*R*W3*dR + R*dW1 - dW4*R + dR*W1 - W4*dR - dW2)

dθdz = -im/2*k*(dW[2] + dW[4]*R + W[4]*dRdθ -
        (dRdθ*W[1] + R*dW[1]) -
        (dRdθ*W[3]*R + R*dW[3]*R + R*W[3]*dRdθ))
