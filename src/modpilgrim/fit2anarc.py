'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.3
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
| Last Update:  2021/04/20 (Y/M/D) |
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
    #return np.log(A) + np.log(tfrac1**n) - B*tfrac2
    return np.log(A) + n*np.log(tfrac1) - B*tfrac2
#--------------------#
def log_anarc5(T,A,B,n,T0,Tr=TR_DEFAULT):
    if A  < 0.0: return np.inf
    if Tr < 0.0: return np.inf
    if T0 < 0.0: return np.inf
    tfrac1 = (T+T0)/Tr
    tfrac2 = (T+T0)/(T**2 + T0**2)
    #return np.log(A*tfrac1**n) - B*tfrac2
    return np.log(A) + n*np.log(tfrac1) - B*tfrac2
    #return np.log(A) + np.log(tfrac1**n) - B*tfrac2
#--------------------#
#def log_anarc1(T,A,B)                   : return np.log(anarc1(T,A,B))
##--------------------#
#def log_anarc2(T,A,B,n)                 : return np.log(anarc2(T,A,B,n))
##--------------------#
#def log_anarc3(T,A,B,n,Tr=TR_DEFAULT)   : return np.log(anarc3(T,A,B,n,Tr))
##--------------------#
#def log_anarc4(T,A,B,n,T0,Tr=TR_DEFAULT): return np.log(anarc4(T,A,B,n,T0,Tr))
##--------------------#
#def log_anarc5(T,A,B,n,T0,Tr=TR_DEFAULT): return np.log(anarc5(T,A,B,n,T0,Tr))
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
# value: (function,number of parameters)
ANFNC    = {}
ANFNC[1] = [anarc1,2] # Arrhenius
ANFNC[2] = [anarc2,3] # van't Hoff type 1
ANFNC[3] = [anarc3,3] # van't Hoff type 2
ANFNC[4] = [anarc4,4] # Truhlar
ANFNC[5] = [anarc5,4] # Truhlar, version 2
#====================#
LOGANFNC    = {}
LOGANFNC[1] = [log_anarc1,2] # Arrhenius
LOGANFNC[2] = [log_anarc2,3] # van't Hoff type 1
LOGANFNC[3] = [log_anarc3,3] # van't Hoff type 2
LOGANFNC[4] = [log_anarc4,4] # Truhlar
LOGANFNC[5] = [log_anarc5,4] # Truhlar, version 2
#====================#


#====================#
def fit2anarc(tlist,rclist,log=True):
    # initialize dictionary
    dfit  = {}
    # number of points
    npts  = len(tlist)
    if npts != len(rclist): return dfit
    #---------#
    # Guesses #
    #---------#
    # Get guesses for each type
    if npts < 2: return dfit
    x1 = 1.0/tlist[ 0]
    x2 = 1.0/tlist[-1]
    y1 = np.log(rclist[ 0])
    y2 = np.log(rclist[-1])
    B   = (y1-y2)/(x2-x1)
    lnA = y1+B*x1
    A   = exp128(lnA)
    A,B = float(A),float(B) # curve fit DOES NOT admit float128
    # dictionary of guesses for each type of rate constant
    GUESSES    = {}
    GUESSES[1] = [A,B]
    GUESSES[2] = [A,B,0.0]
    GUESSES[3] = [A,B,0.0]
    GUESSES[4] = [A,B,0.0,100.0] # A, B, n, T0
    GUESSES[5] = [A,B,0.0,100.0] # A, B, n, T0
    # log is True or False??
    if log:
       FNCS   = LOGANFNC
       rclist = [float(np.log(k)) for k in rclist]
    else:
       FNCS   = ANFNC
       rclist = [float(ki) for ki in rclist]
    #-----------------------#
    # Fitting and r^2       #
    #-----------------------#
    for atype in sorted(FNCS.keys()):
        function, nparams = FNCS[atype]
        guess = GUESSES[atype]
        if npts < nparams: continue
        # fitting
        try   : popt, pcov  = curve_fit(function, tlist, rclist, p0=guess, maxfev=10000)
        except: continue
        # calculate r^2
        if len(popt) != nparams: continue
        r2 = get_r2(tlist,rclist,function,popt)
        # add Tr to return
        if atype in [3,4,5]: popt = list(popt)+[TR_DEFAULT]
        # save data
        dfit[atype] = (popt,r2)
        # save analytic3 for analytic4 guesses
        if atype == 3: GUESSES[4][0:3] = popt[0:3]
        if atype == 4: GUESSES[5][0:4] = popt[0:4]
    return dfit
#====================#


def crc16():
    data = '''
 1000.0    5.266E+06
 1050.0    1.048E+07
 1100.0    1.959E+07
 1150.0    3.467E+07
 1200.0    5.853E+07
 1250.0    9.474E+07
 1300.0    1.478E+08
 1350.0    2.231E+08
 1400.0    3.269E+08
 1450.0    4.666E+08
 1500.0    6.505E+08
 1550.0    8.875E+08
 1600.0    1.188E+09
 1650.0    1.562E+09
 1700.0    2.020E+09
 1750.0    2.576E+09
 1800.0    3.240E+09
 1850.0    4.024E+09
 1900.0    4.943E+09
 1950.0    6.007E+09
 2000.0    7.230E+09'''
    ltemp, lk = [], []
    for line in data.split("\n"):
        line = line.strip()
        if line == "": continue
        T,k = line.split()
        ltemp.append(float(T))
        lk.append(float(k))
    dfit = fit2anarc(ltemp,lk,log=True)
    lines = []
    for ana in sorted(dfit.keys()):
        popt,r2 = dfit[ana]
        coefs = " ".join(["%9.2E"%coef for coef in popt])
        line  =  "analytic%i  %s"%(ana,coefs)
        r2    = "(1-r^2 = %5.0E)"%(1.0-r2)
        lines.append((line,r2))
    ml = max([len(line) for line,r2 in lines])
    for line,r2 in lines:
        print(("k(CRC16.fw)  %%-%is  %%s"%ml)%(line,r2))

if __name__ == '__main__':
   def example(log=False):
       tlist = [1000.00,1100.00,1200.00,1300.00,1400.00,\
                1500.00,1600.00,1700.00,1800.00,1900.00,2000.00 ]
       rclist = [5.578E+09,1.525E+10,3.532E+10,7.198E+10,\
                 1.326E+11,2.253E+11,3.584E+11,5.401E+11,\
                 7.777E+11,1.078E+12,1.446E+12]
       dfit = fit2anarc(tlist,rclist,log)
       for ana in sorted(dfit.keys()):
           popt,r2 = dfit[ana]
           print("analytic%i   "%ana + " ".join(["%9.2E"%coef for coef in popt]))

   print("CRC16:")
   crc16()
   print("")

   print("True:")
   example(True)
   print("")

   print("False:")
   example(False)
   print("")





