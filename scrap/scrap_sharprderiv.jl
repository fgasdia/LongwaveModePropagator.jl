using SymPy

#==
tmatrix
==#

@vars θ M11 M12 M13 M21 M22 M23 M31 M32 M33
S = sin(θ)
C = cos(θ)
C² = C^2

# Denominator of most of the entries of `T`
den = inv(1 + M33)

M31den = M31*den
M32den = M32*den

T11 = -S*M31den
T12 = S*M32den
# T13 = 0
T14 = (C² + M33)*den
# T21 = 0
# T22 = 0
# T23 = 1
# T24 = 0
T31 = M23*M31den - M21
T32 = C² + M22 - M23*M32den
# T33 = 0
T34 = S*M23*den
T41 = 1 + M11 - M13*M31den
T42 = M13*M32den - M12
# T43 = 0
T44 = -S*M13*den

dT11 = diff(T11, θ)
dT12 = diff(T12, θ)
dT14 = diff(T14, θ)
dT31 = diff(T31, θ)
dT32 = diff(T32, θ)
dT34 = diff(T34, θ)
dT41 = diff(T41, θ)
dT42 = diff(T42, θ)
dT44 = diff(T44, θ)

#==
Field coeffs
==#

# We work backwards. dR first
u11 = sympy.Function("u11")(θ)
u21 = sympy.Function("u21")(θ)
u12 = sympy.Function("u12")(θ)
u22 = sympy.Function("u22")(θ)
d11 = sympy.Function("d11")(θ)
d21 = sympy.Function("d21")(θ)
d12 = sympy.Function("d12")(θ)
d22 = sympy.Function("d22")(θ)

d = Sym[d11 d12; d21 d22]
u = Sym[u11 u12; u21 u22]

R = d*(u^-1)

simplify(diff(R[1,1], θ))

# What are du and dd?
s1 = sympy.Function("s1")(θ)
s2 = sympy.Function("s2")(θ)
s3 = sympy.Function("s3")(θ)
s4 = sympy.Function("s4")(θ)

d11 = C*s4-s1
d21 = -C*s2+s3
# d12 = C*s42-s12
# d22 = -C*s22+s32

u11 = C*s4+s1
u21 = -C*s2-s3
# u21 = C*s42+s12
# u22 = -C*s22-s32

dd11 = diff(d11, θ)
dd21 = diff(d21, θ)

du11 = diff(u11, θ)
du21 = diff(u21, θ)

# What are the s's?
q = sympy.Function("q")(θ)

a1 = sympy.Function("a1")(θ)
a2 = sympy.Function("a2")(θ)
a3 = sympy.Function("a3")(θ)
a4 = sympy.Function("a4")(θ)
a5 = sympy.Function("a5")(θ)
a6 = sympy.Function("a6")(θ)
A = sympy.Function("A")(θ)

# Using s instead of e
s1 = a3*q + a4
s2 = A
s3 = q*A
s4 = a5*q + a6

ds1 = diff(s1, θ)
ds2 = diff(s2, θ)
ds3 = diff(s3, θ)
ds4 = diff(s4, θ)

# What is A?
A = q^2 + a1*q + a2
dA = diff(A, θ)

# What are the as?
T11 = sympy.Function("T11")(θ)
T12 = sympy.Function("T12")(θ)
T14 = sympy.Function("T14")(θ)
T41 = sympy.Function("T41")(θ)
T42 = sympy.Function("T42")(θ)
T44 = sympy.Function("T44")(θ)

a1 = -(T11 + T44)
a2 = T11*T44 - T14*T41
a3 = T12
a4 = T14*T42 - T12*T44
a5 = T42
a6 = T12*T41 - T11*T42

da1 = diff(a1, θ)
da2 = diff(a2, θ)
da3 = diff(a3, θ)
da4 = diff(a4, θ)
da5 = diff(a5, θ)
da6 = diff(a6, θ)



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


#==
W matrix deriv
==#
@vars θ

T11 = sympy.Function("T11")(θ)
T12 = sympy.Function("T12")(θ)
T14 = sympy.Function("T14")(θ)
T31 = sympy.Function("T31")(θ)
T32 = sympy.Function("T32")(θ)
T34 = sympy.Function("T34")(θ)
T41 = sympy.Function("T41")(θ)
T42 = sympy.Function("T42")(θ)
T44 = sympy.Function("T44")(θ)

C = cos(θ)
Cinv = 1/cos(θ)


# Precompute
T12Cinv = T12*Cinv
T14Cinv = T14*Cinv
T32Cinv = T32*Cinv
T34Cinv = T34*Cinv
CT41 = C*T41

#==
W = | W11 | W12 |
    | W21 | W22 |

W11 = | a11+a11r | -b11 |
      | -c11     | d11  |

W12 = | a21+a21r | -b11 |
      | c12      | d12  |

W21 = | a21-a21r | b21  |
      | c11      | -d12 |

W22 = | a11-a11r | b21  |
      | -c12     | -d11 |
==#

a11 = T11 + T44
a11r = T14Cinv + CT41
a21 = T44 - T11
a21r = T14Cinv - CT41

b11 = T12Cinv + T42
b21 = T12Cinv - T42

c11 = T31 + T34Cinv
c12 = T31 - T34Cinv

d11 = C + T32Cinv
d12 = T32Cinv - C

# Form the four 2x2 submatrices of `S`
W11 = SMatrix{2,2}(a11+a11r, -c11, -b11, d11)
W12 = SMatrix{2,2}(a21+a21r, c12, -b11, d12)
W21 = SMatrix{2,2}(a21-a21r, c11, b21, -d12)
W22 = SMatrix{2,2}(a11-a11r, -c12, b21, -d11)

dW11a = diff(a11+a11r, θ)
dW11b = diff(-c11, θ)
dW11c = diff(-b11, θ)
dW11d = diff(d11, θ)

dW12a = diff(a21+a21r, θ)
dW12b = diff(c12, θ)
dW12c = diff(-b11, θ)
dW12d = diff(d12, θ)

dW21a = diff(a21-a21r, θ)
dW21b = diff(c11, θ)
dW21c = diff(b21, θ)
dW21d = diff(-d12, θ)

dW22a = diff(a11-a11r, θ)
dW22b = diff(-c12, θ)
dW22c = diff(b21, θ)
dW22d = diff(-d11, θ)
