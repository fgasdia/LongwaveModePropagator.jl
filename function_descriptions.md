# ModeFinder Function Descriptions

## MF_DRIVER
"""
MF_DRIVER sets omega, wavenr, gmax, tmesh, lub, ranger, rangei, and calls MF_WVGD.

LWPC's MF_DRIVER is a short subroutine called by LWP_DRIVER to get starting solutions for modesearch in an individual path segment. I believe _starting solutions_ are the approximate eigenangles. If eigenangles (_modes_) are found, then LWP_DRIVER proceeds to enter SW_STEP and additional routines, presumably to find the exact or complete solutions. 
"""

## MF_WVGD
"""
LWPC's MF_WVGD is a medium length subroutine that obtains approximate eigenangles (_modes_) in a single path segment. It calls several secondary routines to accomplish this task.
"""

## MF_INITEV
"""
MF_INITEV initializes the eigenangle search by calculating the profile, some basic magnetoionic parameters, and calling MF_BRNCPT if necessary. It also determines htntrp.

LWPC's MF_INITEV is a short subroutine. It determines `htntrp`, which is sometimes used as a reference height. It evaluates the profile and calcultes `x`, `y`, `z`, and `g`. If `g` is not low, it calls MF_BRNCPT to locate any branch points in the search grid.

From NOSC TR1143:
- `x = cx*en[1]` is the normalized electron density
- `y = -cy*bfield/omega` is the normalized magnetic field strength
- `z = nu[1]/omega` is ``\nu/\omega`` (collision frequency/angular wave frequency)
- `g` is ``y/(z+im)`` (the ratio of magnetic field to collision frequency)

Some additional values used in MF_INITEV that are not obvious:
- `cx` is in Morfitt Shellman 1976 pg 77 eq 11, which references Wait's NBS TN 300 where
```math
\omega_r(z) = 2.5e5 exp(\beta(z-h\prime)) = 3.182357e9 N(z)/\nu(z)
```
and `cx=3.182357e9`:
```math
\omega_0^2 = 3.18e9 N(0)
\omega = (4\pi N e^2 / m)^(1/2)
```
``\omega_0`` is the plasma frequency, ``N(z)`` is the electron density in electrons/cubic centimeter, and ``\nu(z)`` is the collision frequency in collisions/second.
- `cy = 1.758796e7` is derived (by me) from the ``Y`` in Morfitt Shellman 1976 pg 14 where ``Y = \mu_0\omega_m/\omega`` henries where ``\omega_m`` is the magnetic gyrofrequency (1/s) so that `cy` here is
```math
-cy*bfield/omega = \mu_0 \omega_m / omega = \mu_0 (eB/mc) / omega
cy = \mu_0 (e/mc) * (B/omega)
```
except that 1.758796e7 is just e/mc, so it appears to have no \mu_0.
"""

## MF_BRNCPT
"""
Finds branch points of \zeta and appends them to a list.

At branch points in the theta plane, the ordinary wave cannot be distinguished from the extraordinary wave and magnetoionic reflection coefficients are not defined. This subroutine is a closed solution that calls no other subroutines. It is described in detail in section IV of NOSC TR1143.
"""

## MF_BNDRYS
"""
Builds search grid, including increased resolution around possible Brewster modes.

Medium length subroutine that is self contained math. Not exactly sure what this does. It doesn't exist exactly as-is in Morfitt Shellman 1976, but pg 81 does describe the general process (wihout mention of Brewster modes). It is not clear to me whether or not these boundaries are rectangular or triangular meshes.
"""

## MF_WHICH
"""
Determine if the mesh boundary does **not** contain an eigenangle?

Short subroutine which is entirely logic (no math) and standalone. Seems to isolate a few special cases in which the `noevct` flag should be set, which likely suggests there are no eigenangles here.
"""

## MF_INITOE *important
"""
Determines coefficients for interpolation series of elements of magnetoionic reflection matrix.

Long subroutine that contains math, logic, and calls several other routines including MF_INTEG, MF_EMTRX, MF_FSINTG, MF_ROEMTX, and MF_RMTRX. Details are in section V of NOSC TR 1143, but the subroutine itself is already well commented.

TODO: This is a candidate for replacement with a standard function, or at least partially replaced. 
"""

---

