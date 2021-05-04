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
| Sub-module :  cvtst              |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains functions for
the calculation of the TST rate constant
and the CVT correction factor
'''

#===============================================#
import numpy     as np
#-----------------------------------------------#
import common.fncs       as fncs
import common.Exceptions as Exc
#-----------------------------------------------#
import modpilgrim.steepdesc as sd
#-----------------------------------------------#
from   common.physcons  import KB, KCALMOL, VOL0
from   common.Spline    import Spline
#===============================================#


#===============================================#
def get_gibbs_matrix(dMols,points,VadiSpl,temps,case="cc"):
    '''
    dMOls, a dictionary
        key: label of the point
        value: a tuple (s_i,Mol_i), where s_i is the
               path coordinate and Mol_i is a Molecule
               instance
    '''
    # Initialize matrix with values
    nrows = len(points)
    ncols = len(temps)
    mgibbs = np.zeros( (nrows,ncols) )
    # Calculation of partition functions along MEP
    gibbs   = []
    gibbsTS = None
    for row,point in enumerate(points):
        s_i, Mol_i = dMols[point]
        qtot, V1, (phtra,pfrot,pfvib,pfele) = Mol_i.calc_pfns(temps,case=case,fmode=-1)
        # calculate gibbs free energy (without sigma_rot)
       #rotsigma = Mol_i._rotsigma
       #gibbs_i  = [V1-KB*T*np.log(qtot[idx]*rotsigma) for  idx,T in enumerate(temps)]
        gibbs_i  = [V1-KB*T*np.log(qtot[idx]*VOL0) for  idx,T in enumerate(temps)]
        # append data
        gibbs.append( gibbs_i )
        # save values for TS
        if s_i == 0.0: gibbsTS = list(gibbs_i)
    # add data to mgibbs
    for row in range(nrows):
        for col in range(ncols):
            mgibbs[row][col] = float(gibbs[row][col])
            if gibbsTS is not None: mgibbs[row][col] -= gibbsTS[col]
    # return data
    return mgibbs, gibbsTS
#-----------------------------------------------#
def relative_gamma(xvals,yvals,T):
    # better values in kcal/mol
    yvals = [ yi*KCALMOL for yi in yvals]
    # find maximum
    spl = Spline(xvals,yvals,tan="findiff",tension=0.0)
    spl.find_xtr("max")
    CVT_s, CVT_gibbs = spl.get_max()
    # undo kcal/mol
    CVT_gibbs = CVT_gibbs / KCALMOL
    # CVT values
    CVT_s     = float(CVT_s)
    CVT_gamma = fncs.exp128( -(CVT_gibbs/KB/T))
    # Correct value, just in case
    if CVT_gamma > 1.0: CVT_gamma = 1.0
    # data with more points (for plotting)
    npts = 20*len(xvals)
    dx = (xvals[-1]-xvals[0])/(npts)
    new_xx = [xvals[0]+ii*dx for ii in range(npts+1)]
    new_yy = [spl(x)/KCALMOL for x in new_xx]
    # return
    return float(CVT_s), float(CVT_gamma), (new_xx,new_yy)
#-----------------------------------------------#
def get_cvt(dMols,points,VadiSpl,temps,useics="no"):

    if useics in ["yes",True]: case = "ic"
    else                     : case = "cc"
    xvals = [dMols[point][0] for point in points]
    lcvt_s, lcvt_gamma = [], []

    # get gibbs matrix
    mgibbs, gibbsTS = get_gibbs_matrix(dMols,points,VadiSpl,temps,case)

    # For each temperature, find maximum gibbs
    lnew = []
    for col,T in enumerate(temps):
        yvals = mgibbs[:,col]
        CVT_s, CVT_gamma, new = relative_gamma(xvals,yvals,T)
        # Save data
        lcvt_s.append(CVT_s)
        lcvt_gamma.append(CVT_gamma)
        lnew.append(new)

    # Return data
    return lcvt_s, lcvt_gamma, mgibbs, gibbsTS, lnew
#===============================================#



