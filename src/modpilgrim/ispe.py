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
| Sub-module :  ispe               |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#==============================================#
import numpy     as     np
#----------------------------------------------#
import common.fncs       as  fncs
import common.physcons   as  pc
import common.files      as  ff
import common.Exceptions as  Exc
#----------------------------------------------#
import modpilgrim.steepdesc as     sd
#----------------------------------------------#
from   common.Molecule  import Molecule
from   common.Spline    import Spline
from   common.criteria  import EPS_DLEVELS
#----------------------------------------------#
from   modpilgrim.steepdesc import TSLABEL
#==============================================#

def read_ispe(filename):
    '''
    read input file for ispe
    '''
    # read lines
    lines = ff.read_file(filename)
    lines = fncs.clean_lines(lines,strip=True)
    # initialize data
    ispe_xy = []
    VR , VP = None, None
    tension = 0.0
    # find data in lines
    for line in lines:
        if line == "\n": continue
        label, val = line.split()
        val = float(val)
        if   label.lower() == "tension": tension = val
        elif label.lower() == "reac"   : VR = val
        elif label.lower() == "prod"   : VP = val
        else                           : ispe_xy.append( (label,val) )
    return ispe_xy, tension, VR, VP

#==============================================#
def get_s0_L(VR,VTS,VP,mu,ifreq):
    '''
    calculates s0 and L in base of
    low-level parameters
    '''
    sA0 = -np.sqrt( (VTS-VR)/(ifreq**2)/mu )
    sB0 = +np.sqrt( (VTS-VP)/(ifreq**2)/mu )
    
    sA  = -min(  abs(sA0),2*sB0)
    sB  = +min(2*abs(sA0),  sB0)

    s0 = (sA+sB)/2.0
    L  = (abs(sA)+sB)/2.0
    return s0, L
#----------------------------------------------#
def s2z(sval,s0,L):
    '''
    converts s --> z 
    '''
    arg = (sval-s0)/L
    z  = 2.0/pc.PI * np.arctan(arg)
    return z
#----------------------------------------------#
def get_ts_ifreq(xcc,gcc,Fcc,E,tcommon):
    ch, mtp, atonums, masses, mu = tcommon
    ts = Molecule()
    ts.setvar(xcc=xcc,gcc=gcc,Fcc=Fcc)
    ts.setvar(atonums=atonums,masses=masses)
    ts.setvar(ch=ch,mtp=mtp,V0=E)
    ts.prepare()
    ts.setup(mu)
    ts.ana_freqs()
    [(ifreq, evec)] = ts._ccimag
    return ifreq, mu
#----------------------------------------------#
def in_interval(si,s0,sn):
    return s0-EPS_DLEVELS <= si <= sn+EPS_DLEVELS
#----------------------------------------------#
def gen_data_for_spline(ispe_xy,drst,s0,L):
    lz, ldE = [], []
    for sval,point,Vhl in ispe_xy:
        # low-level data
        Vll = drst[point][1]
        # convert to z-val
        z = s2z(sval,s0,L)
        # get difference
        dE = Vll-Vhl
        # append data
        lz.append(z)
        ldE.append(dE)
    return lz, ldE
#==============================================#

def ispe(tcommon,drst,ispe_xy,tension):
    '''
    V^LL --> V^HL
    '''

    # points in rst
    points = sd.sorted_points(drst,hess=False)

    # Get imaginary frequency for TS (low-level)
    sTS, ETS_ll, xcc, gcc, Fcc = drst[sd.TSLABEL][0:5]
    ifreq, mu = get_ts_ifreq(xcc,gcc,Fcc,ETS_ll,tcommon)
    
    # convert x in ispe_xy to labels
    for idx, (l_i,V_i) in enumerate(ispe_xy):
        if l_i in points:
           s_i = drst[l_i][0]
           label_i = l_i
        else:
           s_i = float(l_i)
           label_i = None
           for point in points:
               if abs(drst[point][0]-s_i) < EPS_DLEVELS: label_i = point
           if label_i is None: raise Exc.DLEVELsthWrong(Exception)
        ispe_xy[idx] = (s_i,label_i,V_i)
    # sort points
    ispe_xy.sort()

    # only one point
    if len(ispe_xy) == 1:
       sx, lx, Ex_hl = ispe_xy[0]
       Ex_ll = drst[lx][1]
       # calculate energy difference
       diffE = Ex_ll - Ex_hl
       # points
       points = sd.sorted_points(drst,hess=False)
       # save old energies
       xx   = [ drst[point][0] for point in points]
       yyll = [ drst[point][1] for point in points]
       # apply difference to all points
       for point in points:
           drst[point] = list(drst[point])
           drst[point][1] = drst[point][1] - diffE
           drst[point] = tuple(drst[point])
    # more than one point
    else:
       # reduce points of drst to those for DLEVEL
       s1, l1, E1_hl = ispe_xy[0]
       sn, ln, En_hl = ispe_xy[-1]
       drst   = {point:drst[point] for point in points if in_interval(drst[point][0],s1,sn)}
       points = sd.sorted_points(drst,hess=False)

       # get s0 and L from lower level
       E1_ll, En_ll = drst[l1][1], drst[ln][1]
       s0, L = get_s0_L(E1_ll,ETS_ll,En_ll,mu,ifreq)

       # generate data for spline
       lxx, lyy = gen_data_for_spline(ispe_xy,drst,s0,L)

       # create spline
       spl = Spline(lxx,lyy,tension=tension)

       # save old energies
       xx   = [ drst[point][0] for point in points]
       yyll = [ drst[point][1] for point in points]

       # modify drst
       for point in points:
           s, Vll =  drst[point][0:2]
           z      = s2z(s,s0,L)
           drst[point] = list(drst[point])
           drst[point][1] = Vll - spl(z)
           drst[point] = tuple(drst[point])

    # save new energies
    yyhl = [ drst[point][1] for point in points]

    # return data
    return drst, points, xx, yyll, yyhl

