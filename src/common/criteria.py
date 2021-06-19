'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.4
License     : MIT/x11

Copyright (c) 2021, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Module     :  common             |
| Sub-module :  criteria           |
| Last Update:  2020/04/18 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''


# |sin(x)|   < EPS_SCX --> sin(x) =  0
# |cos(x)-1| < EPS_SCX --> cos(x) =  1
# |cos(x)+1| < EPS_SCX --> cos(x) = -1
EPS_SCX = 1e-7

# small angle
EPS_SMALLANGLE = 1e-3

# scale for the determination of molecular connectivity
CONNECTSCAL = 1.3

# angle in (180-EPS_LINEAR,180+EPS_LINEAR) is considered linear
EPS_LINEAR = 4.5

# ds STEP FOR THE CALCULATION OF RETURN POINTS
DS_RPT = 1e-3

# comparision of geometries (bohr)
EPS_GEOM = 1e-5

# zero moment of inertia
EPS_INERTIA = 1e-7

# the norm of a vector is ZERO
EPS_NORM = 1e-7

# step for numerical calculation of hessian matrix in analitic surfaces
EPS_HESSDX = 1e-4

# step for cubic first step in MEP calculation 
STEP_CUBIC = 1e-4

# comparision of float numbers
EPS_FLOAT = 1e-8

# comparision of float numbers (strict)
EPS_FSTRICT = 1e-10

# comparision of temperatures
EPS_TEMP = 1e-8

# Small energy difference in kcalmol
EPS_KCALMOL = 0.01

# zero por single-value decomposition
EPS_SVD = 1e-9

# zero por generalized inverse
EPS_GIV = 1e-10

# Ignore freqs smaller than EPS_IC cm^-1 when cleaning frequencies
EPS_IC = 0.3

# Ignore ic-freqs smaller than EPS_ICF cm^-1
EPS_ICF = 0.3

# If cc-freq and ic-freq differ less than EPS_CCIC cm^-1, they can be considered to be equal
EPS_CCIC = 0.5

# In kcal/mol
EPS_CCICZPE = 0.1

# comparison of s values of MEP (in bohr)
EPS_MEPS = 1e-7

# comparison of energies of MEP (in hartree)
EPS_MEPE = 1e-8

# Criteria to increase MEP for one or both sides (in kcal/mol)
EPS_MEPINCR = 1.5

# to identify s values in DLEVEL
EPS_DLEVELS = 1e-9

# zero for the comparison of masses in amu
EPS_AMU = 0.001

