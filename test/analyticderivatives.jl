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


@vars q B4 B3 B2 B1 B0
