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
| Module     :  modpilgrim         |
| Sub-module :  fit2anarc          |
| Last Update:  2021/11/22 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different functions
to adjust data (temperatures and rate
constants) to analytic functions.
'''

#=======================================================#
import numpy as np
from   common.physcons import KB
from   common.fncs     import exp128
from   scipy.optimize  import curve_fit
#-------------------------------------------------------#
np.warnings.filterwarnings('ignore')
#=======================================================#

TR_DEFAULT = 300.0

#====================#
# Analytic functions #
#====================#
def anarc1(T,A,B):
    '''
    Arrhenius: k = A*exp(-B/T);
      * A in cm^3/molecule/s
      * B in K.
    '''
    # to avoid "RuntimeWarning: invalid value encountered in power"
    if A  < 0.0: return np.inf
    # Calculate rate constant
    k = A * exp128(-B/T)
    # Return k
    return k
#--------------------#
def anarc2(T,A,B,n):
    '''
    Van't Hoff type 1: k = A*T^n*exp(-B/T)
      * A*T^n in cm^3/molecule/s.
      * B in K.
      * n (adimensional)
    '''
    # to avoid "RuntimeWarning: invalid value encountered in power"
    if A  < 0.0: return np.inf
    # Calculate rate constant
    k = A * (T**n) * exp128(-B/T)
    # Return k
    return k
#--------------------#
def anarc3(T,A,B,n,Tr=TR_DEFAULT):
    '''
    Van't Hoff type 2: k = A*(T/Tr)^n*exp(-B/T)
      * A in cm^3/molecule/s.
      * B in K.
      * n (adimensional)
      * Tr in Kelvin
    '''
    # to avoid "RuntimeWarning: invalid value encountered in power"
    if A  < 0.0: return np.inf
    if Tr < 0.0: return np.inf
    # Calculate rate constant
    k = A * ((T/Tr)**n) * exp128(-B/T)
    # Return k
    return k
#--------------------#
def anarc4(T,A,B,n,T0,Tr=TR_DEFAULT):
    '''
    Truhlar: k = A*(T/Tr)^n*exp(-B*(T+T0)/(T^2+T0^2))
      * A in cm^3/molecule/s.
      * B in K.
      * n (adimensional)
      * Tr in Kelvin
      * T0 in Kelvin
    '''
    # to avoid "RuntimeWarning: invalid value encountered in power"
    if A  < 0.0: return np.inf
    if Tr < 0.0: return np.inf
    # avoid negative values
    if T0 < 0.0: return np.inf
    # Calculate rate constant
    tfrac1 = T/Tr
    tfrac2 = (T+T0)/(T**2 + T0**2)
    k = A * (tfrac1**n) * exp128(-B*tfrac2)
    # Return k
    return k
#--------------------#
def anarc5(T,A,B,n,T0,Tr=TR_DEFAULT):
    '''
    Truhlar: k = A*((T+T0)/Tr)^n*exp(-B*(T+T0)/(T^2+T0^2))
      * A in cm^3/molecule/s.
      * B in K.
      * n (adimensional)
      * Tr in Kelvin
      * T0 in Kelvin
    '''
    # to avoid "RuntimeWarning: invalid value encountered in power"
    if A  < 0.0: return np.inf
    if Tr < 0.0: return np.inf
    # avoid negative values
    if T0 < 0.0: return np.inf
    # Calculate rate constant
    tfrac1 = (T+T0)/Tr
    tfrac2 = (T+T0)/(T**2 + T0**2)
    k = A * (tfrac1**n) * exp128(-B*tfrac2)
    # Return k
    return k
#====================#

#====================#
# Logarithm versions #
#====================#
def log_anarc1(T,A,B):
    if A  < 0.0: return np.inf
    return np.log(A) - B / T
#--------------------#
def log_anarc2(T,A,B,n):
    if A  < 0.0: return np.inf
    return np.log(A * (T**n)) - B/T
#--------------------#
def log_anarc3(T,A,B,n,Tr=TR_DEFAULT):
    if A  < 0.0: return np.inf
    if Tr < 0.0: return np.inf
    return np.log(A * ((T/Tr)**n)) - B/T
#--------------------#
def log_anarc4(T,A,B,n,T0,Tr=TR_DEFAULT):
    if A  < 0.0: return np.inf
    if Tr < 0.0: return np.inf
    if T0 < 0.0: return np.inf
    tfrac1 = T/Tr
    tfrac2 = (T+T0)/(T**2 + T0**2)
    return np.log(A) + n*np.log(tfrac1) - B*tfrac2
#--------------------#
def log_anarc5(T,A,B,n,T0,Tr=TR_DEFAULT):
    if A  < 0.0: return np.inf
    if Tr < 0.0: return np.inf
    if T0 < 0.0: return np.inf
    tfrac1 = (T+T0)/Tr
    tfrac2 = (T+T0)/(T**2 + T0**2)
    return np.log(A) + n*np.log(tfrac1) - B*tfrac2
#====================#


#=====================#
# Activation Energies #
#=====================#
def activation1(T,A,B):
    return B*KB
#--------------------#
def activation2(T,A,B,n):
    return (B+n*T) * KB
#--------------------#
def activation3(T,A,B,n,Tr=TR_DEFAULT):
    return (B+n*T) * KB
#--------------------#
def activation4(T,A,B,n,T0,Tr=TR_DEFAULT):
    temp_term = (T**4 + 2*T0*(T**3) - (T0**2)*(T**2)) / ((T**2 + T0**2)**2)
    return (B*temp_term+n*T) * KB
#--------------------#
def activation5(T,A,B,n,T0,Tr=TR_DEFAULT):
    temp_term = (T**4 + 2*T0*(T**3) - (T0**2)*(T**2)) / ((T**2 + T0**2)**2)
    return (B*temp_term + n*(T*T/(T+T0))) * KB
#=====================#



#====================#
# Get r^2 of fitting #
#====================#
def get_r2(xdata,ydata,function,params):
    # mean value of y
    ymean = sum(ydata)/len(ydata)
    # fitted values
    yfit  = [function(xi,*params) for xi in xdata]
    # get r^2
    SStot,SSres = 0.0, 0.0
    for y,yf in zip(ydata,yfit):
        SSres += (y-yf)**2
        SStot += (y-ymean)**2
    rsquare = 1.0- SSres/SStot
    # return r^2
    return rsquare
#====================#


#====================#
# value: (function,ln(function),number of parameters)
#====================#
FNCS    = {}
FNCS[1] = [anarc1,log_anarc1,2] # Arrhenius
FNCS[2] = [anarc2,log_anarc2,3] # van't Hoff type 1
FNCS[3] = [anarc3,log_anarc3,3] # van't Hoff type 2
FNCS[4] = [anarc4,log_anarc4,4] # Truhlar, version 1
FNCS[5] = [anarc5,log_anarc5,4] # Truhlar, version 2
#====================#


#====================#
def fit2anarc(tlist,rclist,log=True):
    # initialize variables
    dfit  = {}         # dictionary
    npts  = len(tlist) # number of points
    # check consistency in number of points
    if npts != len(rclist): return dfit
    #---------#
    # Guesses #
    #---------#
    # Get guesses for each type
    if npts < 2: return dfit
    x1,y1 = 1.0/tlist[ 0], np.log(rclist[ 0])
    x2,y2 = 1.0/tlist[-1], np.log(rclist[-1])
    B     = (y1-y2)/(x2-x1)
    lnA   = y1+B*x1
    A,B   = float(exp128(lnA)),float(B) # curve fit DOES NOT admit float128
    # dictionary of guesses for each type of rate constant
    GUESSES    = {}
    GUESSES[1] = [A,B]
    GUESSES[2] = [A,B,0.0]
    GUESSES[3] = [A,B,0.0]
    GUESSES[4] = [A,B,0.0,100.0] # A, B, n, T0
    GUESSES[5] = [A,B,0.0,100.0] # A, B, n, T0

    # rate constants and their logarithm
    list_k   = [float(       k ) for k in rclist]
    list_lnk = [float(np.log(k)) for k in rclist]

    #-----------------------#
    # Fitting and r^2       #
    #-----------------------#
    for atype in sorted(FNCS.keys()):

        function, lnfunction, nparams = FNCS[atype]
        if npts < nparams: continue

        # what to fit: k or log(k)?
        if log: the_function, the_ydata = lnfunction, list_lnk
        else  : the_function, the_ydata =   function, list_k
        # do fitting
        guess = list(GUESSES[atype])
        try   : popt,pcov = curve_fit(the_function, tlist, the_ydata, p0=guess, maxfev=10000)
        except: continue

        # calculate r^2 considering k (and never log(k))
        if len(popt) != nparams: continue
        r2 = get_r2(tlist,list_k,function,popt)

        # add Tr to optimized parameters
        if atype in [3,4,5]: popt = list(popt)+[TR_DEFAULT]

        # save data
        dfit[atype] = (popt,r2)

        # use parameters as guesses for next fittings
#       if atype == 3: GUESSES[4][0:3] = popt[0:3]
#       if atype == 4: GUESSES[5][0:4] = popt[0:4]
    return dfit
#====================#


