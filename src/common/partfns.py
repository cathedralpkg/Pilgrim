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
| Sub-module :  partfns            |
| Last Update:  2020/06/04 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains some functions related
to the calculation of basic partition functions

'''

#===============================================#
import numpy           as     np
#-----------------------------------------------#
import common.physcons as     pc
from   common.fncs     import exp128
#===============================================#

#===============================================#
# Partition functions                           #
#===============================================#
def pf_partinbox(mass,T):
    return (pc.TWOPI*pc.KB*mass*T)**(3./2.)/(pc.H**3)
#-----------------------------------------------#
def pf_rigidrotor(imoments,T,rotsymnum=1):
    beta      = pc.KB * T
    rot_const = [(pc.HBAR**2)/Ii/2 for Ii in imoments]
    # linear case
    if len(imoments) == 1: qrr = beta / rot_const[0]
    else                 : qrr = np.sqrt(pc.PI) * np.sqrt(beta**3 / np.prod(rot_const))
    return qrr/rotsymnum
#-----------------------------------------------#
def pf_harmosc1D(angfreq,T,imag=1E10):
    if   angfreq  < 0.0: return imag
    elif angfreq == 0.0: return 1.0
    exp  = exp128(-pc.HBAR*angfreq/pc.KB/T)
    qHO  = 1.0/(1.0-exp)
    return qHO
#-----------------------------------------------#
def pf_harmosc(angfreqs,T,imag=1E10):
    qHO = np.prod([pf_harmosc1D(angfreq,T,imag) for angfreq in angfreqs])
    return qHO
#-----------------------------------------------#
def pf_electr(eslist,T):
    pf, beta = 0.0, 1.0 / pc.KB / T
    for deg, relE in eslist: pf += deg * exp128(-relE * beta)
    return pf
#===============================================#

#===============================================#
#      Partition functions and derivatives      #
#===============================================#
def pfnder_tr(mass,T):
    ph_tra   = pf_partinbox(mass,T)
    # at constant volume
    fdln_trV = -1.5*pc.KB*T
    sdln_trV = +1.5*(pc.KB*T)**2
    # at constant pressure
    fdln_trP = -2.5*pc.KB*T
    sdln_trP = +2.5*(pc.KB*T)**2
    return ph_tra,fdln_trV,sdln_trV,fdln_trP,sdln_trP
#-----------------------------------------------#
def pfnder_rot(imoments,T,rotsigma=1):
    pfn_rot  = pf_rigidrotor(imoments,T,rotsigma)
    fdln_rot = -1.5*pc.KB*T
    sdln_rot = +1.5*(pc.KB*T)**2
    return pfn_rot, fdln_rot, sdln_rot
#-----------------------------------------------#
def pfnder_vib_i(freq,T):
    '''
    Calculates vibrational partition function (Q) and
    log(Q) derivatives with regards to beta
    Zero taken as the energy of level n=0 (i.e. the ZPE)
    ------
    T: float or np.array
    '''
    zpe, qvib, fdln, sdln = 0.0, 1.0, 0.0, 0.0 # for imag freqs
    if freq > 0.0:
       zpe    = pc.HBAR * freq / 2.0
       hw     = pc.HBAR * freq
       bhw    = (1.0/pc.KB/T) * hw
       exp    = np.exp(-bhw)
       # partition functions and derivatives (of logQ)
       qvib   = 1.0/(1.0-exp)
       fdln   = - hw * exp / (1.0-exp)
       sdln   = hw**2 / (1/exp - 1.0) / (1.0-exp)
    return zpe, qvib, fdln, sdln
#-----------------------------------------------#
def pfnder_vib(freqs,Tlist):
    zpe      = 0.0
    pfn_vib  = [1.0 for T in Tlist]
    fdln_vib = [0.0 for T in Tlist]
    sdln_vib = [0.0 for T in Tlist]
    for idx,T in enumerate(Tlist):
        for freq in freqs:
            zpe_i, qvib_i, fdln_i, sdln_i = pfnder_vib_i(freq,T)
            # update zpe
            if idx == 0: zpe += zpe_i
            # update partition fncs and derivatives
            pfn_vib[idx]  *= qvib_i
            fdln_vib[idx] += fdln_i
            sdln_vib[idx] += sdln_i
    return zpe, pfn_vib, fdln_vib, sdln_vib
#-----------------------------------------------#
def pfnder_ele(les,Tlist):
    pfn_ele = [0.0 for T in Tlist]
    fd_ele  = [0.0 for T in Tlist]
    sd_ele  = [0.0 for T in Tlist]
    for idx,T in enumerate(Tlist):
        beta = 1.0/pc.KB/T
        for mtp,relE in les:
            pfn_ele[idx] +=           mtp*exp128(-beta*relE)
            fd_ele[idx]  += (relE   )*mtp*exp128(-beta*relE)
            sd_ele[idx]  += (relE**2)*mtp*exp128(-beta*relE)
    # convert to array
    pfn_ele = np.array(pfn_ele)
    fd_ele  = np.array(fd_ele)
    sd_ele  = np.array(sd_ele)
    # we care about derivatives of logQ
    fdln_ele = fd_ele / pfn_ele
    sdln_ele = sd_ele / pfn_ele - (fd_ele/pfn_ele)**2
    return pfn_ele,fdln_ele,sdln_ele
#===============================================#





#==========================================================#
#        Calculation of equilibrium/rate constants         #
#==========================================================#
def Qs2Kc(ltemp,QA,QB,VA,VB,nR=1,nP=1):
    '''
    Qi; partition function per unit volume
    Kc = [B]/[A]
    '''
    term_V = pc.VOL0**(nP-nR)
    return [term_V*QB[idx]/QA[idx]*exp128(-(VB-VA)/pc.KB/T) for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def Kc2GFE(ltemp,Kc):
    return [ -pc.KB * T * np.log(Kc[idx]) for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def Kc2rate(ltemp,Kc):
    return [(pc.KB*T/pc.H)*Kc[idx] for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def rate2Kc(ltemp,k):
    return [k[idx]*pc.H/(pc.KB*T) for idx,T in enumerate(ltemp)]
#----------------------------------------------------------#
def rate2GFE(ltemp,rates,nR=1):
    term_V = pc.VOL0**(1-nR)
    GFE     = [-pc.KB * T * np.log(term_V*ki*pc.H/pc.KB/T) for T,ki in zip(ltemp,rates)]
    return GFE
#==========================================================#


