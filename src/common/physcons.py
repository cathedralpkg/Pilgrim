'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.5
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
| Sub-module :  physcons           |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different physical constants
'''

#==============================================#
# SOME USEFUL CONSTANTS AND CONVERSION FACTORS #
#==============================================#
#------------- related to radians -------------#
PI       = 3.141592653589793                   # from math module
TWOPI    = 6.283185307179586                   #
#---------- physical constants  (SI) ----------#
H_SI     = 6.62606896E-34                      # Planck constant    [J*s]
NA_SI    = 6.022140857E+023                    # Avogradro's number [mol^-1]
C0_SI    = 2.99792458E8                        # speed of light     [m/s]
KB_SI    = 1.3806504E-23                       # Boltzmann const    [J/K]
QE_SI    = 1.6021766208E-019                   # charge of electron [C]
ME_SI    = 9.10938356E-031                     # mass of electron   [kg]
MC12_SI  = 1.66053904E-027                     # mass of C12        [kg]
CAL_SI   = 4.184                               # 1 cal = 4.184 J
R_SI     = KB_SI *  NA_SI                      # Ideal gas constant [J*K^-1*mol^-1]
HBAR_SI  =  H_SI / TWOPI                       # Planck constant divided by 2pi    
#------------- conversion from au -------------#
JOULE    = 4.35974465E-018                     #
METER    = 0.52917721067E-10                   #
SECOND   = 1.0 / 4.1341373337e+16              #
ANGSTROM = METER * 1E10                        #
AMU      = ME_SI/MC12_SI                       #
KG       = ME_SI                               #
CAL      = JOULE / CAL_SI                      #
KCAL     = CAL /1000.0                         #
KJMOL    = JOULE * NA_SI / 1000.0              #
KCALMOL  = KJMOL / CAL_SI                      #
CALMOL   = KCALMOL * 1000.0                    #
JMOL     = KJMOL * 1000.0                      #
METER3   = METER**3                            #
CM       = 100 * METER                         #
ML       = CM**3                               #
EV       = JOULE/QE_SI                         #
DEBYE    = 2.541746                            #
#---------- physical constants in au ----------#
H        = 2*PI                                #
HBAR     = 1.0                                 #
NA       = NA_SI                               #
C0       = C0_SI / (METER/SECOND)              #
KB       = KB_SI / (JOULE)                     #
ME       = ME_SI / (KG)                        #
R        = KB * NA                             #
H2CM     = 1.0 / (H * C0) / CM                 # hartree to 1/cm
CM2H     = 1.0/H2CM                            # 1/cm to hartree
#==============================================#

#==============================================#
#           References states in au            #
#==============================================#
PRE0 = 1.00E5  / (JOULE / METER3)              # reference pressure of 1 bar
VOL0 = 1.00    / ML                            # reference volume of 1 mL per molecule
#==============================================#


