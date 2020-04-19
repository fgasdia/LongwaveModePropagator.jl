using SymPy

########
# Solving eigenvector problem Te = qe => (T- qI)e = 0
# Assuming eigvalues `q` have already been identified, then plug `q` into the
# equation so we're essentially solving Ax = 0
# Below is the analytical solution where we choose e(2) = 1, which means e(3) = q

@vars T11 T12 T14 T31 T32 T34 T41 T42 T44 q E1 E4

# additional equation : T31*e1+T32-q^2+T34*e4
solve([(T11-q)*E1+T12+T14*E4, T41*E1+T42+(T44-q)*E4], (E1, E4))

#==
Dict{Any,Any} with 2 entries:
  e1 => (T12*(T44 - q) - T14*T42)/(T14*T41 - (T11 - q)*(T44 - q))
  e4 => (-T12*T41 + T42*(T11 - q))/(T14*T41 - (T11 - q)*(T44 - q))
==#

# Complete system w/o assuming e(2)=1, e(3)=q
@vars T11 T12 T14 T31 T32 T34 T41 T42 T44 q E1 E2 E3 E4

S = solve([(T11-q)*E1+T12*E2+T14*E4, -q*E1+E3, T31*E1+T32*E2-q*E3+T34*E4, T41*E1+T42*E2+(T44-q)*E4], (E1, E2, E3, E4))
# OK, it doesn't give us a useful answer, only E[:] = 0


# ISOTROPIC

@vars T14 T31 T32 T41 T42 q E1 E2 E3 E4

solve([-q*E1+T14*E4, -q*E2+E3, T31*E1+T32*E2-q*E3, T41*E1+T42*E2-q*E4], (E1, E2, E3, E4))


# Determinant is used to find eigenvalues
D = det([-q 0 0 T14; 0 -q 1 0; T31 T32 -q 0; T41 T42 0 -q])

# Roots q of D
roots(D, q)


########
#==
Sheddy sharp reflection
==#

@vars M11 M12 M13 M21 M22 M23 M31 M32 M33 D11 D12 D13 D21 D22 D23 D31 D32 D33
@vars EX EY EZ q S

L = [-q^2 0 q*S; 0 -q^2-S^2 0; q*S 0 -S^2]
M = [M11 M12 M13; M21 M22 M23; M31 M32 M33]

D = [1 0 0; 0 1 0; 0 0 1] + M + L

Dm = [D11 D12 D13; D21 D22 D23; D31 D32 D33]

E = [EX; EY; EZ]

D*E

# Ey
simplify((D31*D13/(D32*D11) - D33/D32)/(1 - D31*D12/(D32*D11)))

# Ez
simplify(Ez == )
