# level: HF
# Atomic number and non-scaled cartesian coordinates [bohr]
start_cc
   001   +0.00000000E+00  +0.00000000E+00  -1.15549419E-02
   001   +0.00000000E+00  +0.00000000E+00  +1.33436323E+00
end_cc

# Charge, multiplicity, energy [hartree],
# point group and rotational symmetry number
start_basic
   charge        0
   multiplicity  1
   energy       -1.11750590      # Total energy in hartree
   pointgroup    Dinfv           # Point group
   rotsigma      2               # Rotational sigma
end_basic

# Non-scaled cartesian gradient [hartree/bohr]
start_grad
   +0.00000000E+00  +0.00000000E+00  +6.80000000E-07
   +0.00000000E+00  +0.00000000E+00  -6.80000000E-07
end_grad

# Low triangular part of force constant (hessian) matrix [hartree/bohr^2]
# i.e.: F_11, F_21, F_22, F_13, F_23, F_33...
start_hess
   -5.10000000E-07  +0.00000000E+00  -5.10000000E-07  +0.00000000E+00  +0.00000000E+00
   +5.72936510E-01  +5.10000000E-07  +0.00000000E+00  +0.00000000E+00  -5.10000000E-07
   +0.00000000E+00  +5.10000000E-07  +0.00000000E+00  +0.00000000E+00  -5.10000000E-07
   +0.00000000E+00  +0.00000000E+00  -5.72936510E-01  +0.00000000E+00  +0.00000000E+00
   +5.72936510E-01
end_hess

