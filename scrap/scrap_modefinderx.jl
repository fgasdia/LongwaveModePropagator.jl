using SymPy


#==
dTdC
==#
@vars θ x M11 M12 M13 M21 M22 M23 M31 M32 M33
C =  cos(θ)
S = sin(θ)
C² = C^2

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

# d/dcos(θ) (sin(θ)) = -cot(θ)
dT11 = cot(θ)*M31den
dT12 = -cot(θ)*M32den
dT14 = diff(T14.subs(C,x), x).subs(x,C)
dT31 = 0
dT32 = diff(T32.subs(C,x), x).subs(x,C)
dT34 = -cot(θ)*M23*den
dT41 = 0
dT42 = 0
dT44 = cot(θ)*M13*den


#==
dXdCdz
==#

@vars x θ k

T11 = sympy.Function("T11")(C)
T12 = sympy.Function("T12")(C)
T14 = sympy.Function("T14")(C)
T31 = sympy.Function("T31")(C)
T32 = sympy.Function("T32")(C)
T34 = sympy.Function("T34")(C)
T41 = sympy.Function("T41")(C)
T42 = sympy.Function("T42")(C)
T44 = sympy.Function("T44")(C)

X = sympy.Function("X")(C)

C = cos(θ)

a11 = 4*T41
b11 = 2*(T44 - C*T41)
b12 = -2*T42
b22 = -2*C
c11 = 2*(T11 - C*T41)
c21 = 2*T31
c22 = -2*C
d11 = C*(T11 -T44) - T14 + C^2*T41
d21 = -C*T31 + T34
d12 = T12 + C*T42
d22 = C^2 - T32

a = [a11 0; 0 0]
b = [b11 b12; 0 b22]
c = [c11 0; c21 c22]
d = [d11 d12; d21 d22]

da = 0

db11 = 2*(dT44 - T41) = 2*dT44
db12 = 0
db22 = -2

dc11 = 2*(dT11 - T41) = 2*dT11
dc21 = 0
dc22 = -2

# dd11 = diff(d11.subs(C,x), x).subs(x,C)
dd11 = C*(dT11 - dT44) + T11 - T44 + C^2*dT41 + 2*T41*C - dT14
     = C*(dT11 - dT44) + T11 - T44 + 2*T41*C - dT14
dd21 = -C*dT31 - T31 + dT34 = -T31 + dT34
dd12 = dT12 + T42 + C*dT42 = dT12 + T42
dd22 = 2*C - dT32

dX = -1im/2*k*(a + b*X + X*c + X*d*X)

dXdC = -1im/2*k*(da + db*X + b*dX + dX*c + X*dc +
                 dX*d*X + X*dd*X + X*d*dX)

     = -1im/2*k*(db*X + b*dX + dX*c + X*dc + dX*d*X + X*dd*X + X*d*dX)


#==
freespace X wrt C
==#

@vars x θ α k K L n² nb² CH

"""
See MS76 pg 24

zb is the "from" level
"""
function freespaceintegration(z, X, ea::EigenAngle, frequency::Frequency)
    S, C, C² = ea.sinθ, ea.cosθ, ea.cos²θ
    k = frequency.k

    zb = BOTTOMHEIGHT

    α = 2/EARTHRADIUS
    K = 1im*cbrt(α/k)
    L = 1im*(α/(2*K))

    n² = 1 + α*(z - CURVATURE_HEIGHT)
    nb² = 1 + α*(zb - CURVATURE_HEIGHT)

    q = (k/α)^(2/3)*(C² + α*(z - CURVATURE_HEIGHT))
    qb = (k/α)^(2/3)*(C² + α*(zb - CURVATURE_HEIGHT))

    h₁, h₂, h₁p, h₂p = modifiedhankel(q)
    h₁b, h₂b, h₁pb, h₂pb = modifiedhankel(qb)

    a11 = h₁b
    a12 = h₂b
    a21 = C*h₁b + K*h₁pb
    a22 = C*h₂b + K*h₂pb
    Δ₁ = a11*a22 - a21*a12

    AC1 = X[2,1]*a22/Δ₁
    BC1 = -X[2,1]*a21/Δ₁
    AC2 = (X[2,2]*a22 - 2*a12)/Δ₁
    BC2 = (2*a11 - X[2,2]*a21)/Δ₁

    a21 = C*h₁b + (K*h₁pb + L*h₁b)/nb²
    a22 = C*h₂b + (K*h₂pb + L*h₂b)/nb²
    Δ₂ = a11*a22 - a21*a12

    QC1 = (X[1,1]*a22 - 2*a12)/Δ₂
    GC1 = (2*a11 - X[1,1]*a21)/Δ₂
    QC2 = X[1,2]*a22/Δ₂
    GC2 = -X[1,2]*a21/Δ₂

    a21 = C*h₁ + K*h₁p
    a22 = C*h₂ + K*h₂p

    ey1 = (AC1*a21 + BC1*a22)/2
    ey2 = (AC2*a21 + BC2*a22)/2

    a21 = C*h₁ + (K*h₁p + L*h₁)/n²
    a22 = C*h₂ + (K*h₂p + L*h₂)/n²

    ep1 = (QC1*a21 + GC1*a22)/2
    ep2 = (QC2*a21 + GC2*a22)/2

    # Reflection coefficients at `z`
    v1 = AC1*h₁ + BC1*h₂
    v2 = AC2*h₁ + BC2*h₂
    w1 = ep1*ey2 - ep2*ey1

    X21 = (v1*ey2 - v2*ey1)/w1
    X22 = (v2*ep1 - v1*ep2)/w1

    v1 = QC1*h₁ + GC1*h₂
    v2 = QC2*h₁ + GC2*h₂
    w2 = ep1*ey2 - ep2*ey1

    X11 = (v1*ey2 - v2*ey1)/w2
    X12 = (v2*ep1 - v1*ep2)/w2

    return SMatrix{2,2}(X11, X21, X12, X22)
end
