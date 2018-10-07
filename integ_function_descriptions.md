# INTEG steps and functions

## MF_INTEG
"""
Perform an integration of the differential equations for the ionospheric reflection matrix using RK integration formulas. The integration variables are the elements of the matrix (R+1)/C where R is the reflection matrix described by Budden. The (R+1)/C variables used in the comparison step are denoted as `X`.
"""
function MF_INTEG()
# Initialize c, s, csq, ssq

# call MF_SMTRX

end

## MF_SMTRX
"""
Computation of coefficients used in the differential equations for (R+1)/C. Coefficients are derived from the differential equation for the reflection coefficient matrix given by Budden in Proc. Roy. Soc. A227, pp 516-537 (1955).
"""
function MF_SMTRX()
# call MF_TMTX
end

## MF_TMTRX
"""
Computation of M matrix elements as defined by Budden in Proc. Roy. Soc. A227, pp. 516-537 (1955). Thos combinations of M matrix elements used in the T matrix which do not include use of theta are the final output of this routine. The dummy indexing on L is used in view of possible storage (as functions of height) and reuse of these combinations for different values of theta.
"""
function MF_TMTRX()
    # Calculate X, Y, Z, U for each species
    for s = 1:numspecies
        X = N e^2 / \epsilon\_0 m \omega^2
          = N (e/\omega)^2 / m
        Y = \mu_0 e H / m \omega
        Z = \nu / \omega
        U = 1 - im*Z
    end

    M11 = U^2 - dcl^2*Y^2 + curvature
    M22 = U^2 - dcm^2*Y^2 + curvature
    M33 = U^2 - dcn^2*Y^2 + curvature

    M12 = -im*dcn*Y*U - dcl*dcm*Y^2
    M21 = im*dcn*Y*U - dcl*dcm*Y^2

    M13 = im*dcm*Y*U - dcl*dcn*Y^2
    M31 = -im*dcm*Y*U - dcl*dcn*Y^2

    M23 = -im*dcl*Y*U - dcm*dcn*Y^2
    M32 = im*dcl*Y*U - dcm*dcn*Y^2

    M = -X/(U*(U^2-Y^2)) * [M11 12 13; M21 M22 M23; M31 M32 M33]

    S = sin(\theta)
    C = cos(\theta)
    T11 = -S*M31/(1+M33)
    T12 = S*M32/(1+M33)
    T13 = 0
    T14 = (C^2+M33)/(1+M33)
    T21 = 0
    T22 = 0
    T23 = 1
    T24 = 0
    T31 = M23
end
