using SymPy

# Solving eigenvector problem Te = qe => (T- qI)e = 0
# Assuming eigvalues `q` have already been identified, then plug `q` into the
# equation so we're essentially solving Ax = 0
# Below is the analytical solution where we choose e(2) = 1, which means e(3) = q

@vars T11 T12 T14 T31 T32 T34 T41 T42 T44 q e1 e4

# additional equation : T31*e1+T32-q^2+T34*e4
solve([(T11-q)*e1+T12+T14*e4, T41*e1+T42+(T44-q)*e4], (e1, e4))

#==
Dict{Any,Any} with 2 entries:
  e1 => (T12*(T44 - q) - T14*T42)/(T14*T41 - (T11 - q)*(T44 - q))
  e4 => (-T12*T41 + T42*(T11 - q))/(T14*T41 - (T11 - q)*(T44 - q))
==#
