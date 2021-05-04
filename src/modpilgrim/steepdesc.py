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
| Sub-module :  steepdesc          |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different functions
related to the calculation of a steepest
descent path
'''

#===============================================#
import sys
import numpy    as np
import common.fncs     as fncs
import common.internal as intl
import common.Exceptions      as     Exc
from   scipy.integrate import quad
from   scipy.optimize  import brenth
from   common.physcons        import AMU
from   common.criteria        import STEP_CUBIC
from   common.criteria        import EPS_MEPS
from   common.criteria        import EPS_FLOAT
#===============================================#

TSLABEL = "saddle"

#===============================================#
# Basic functions for Steepest descent (sd)     #
#===============================================#
def sorted_points(drst,hess=True):
    if hess:
       points = [(data[0],target) for (target,data) in drst.items()\
                 if data[4] not in [None,[]]]
    else:
       points = [(data[0],target) for (target,data) in drst.items()]
    points.sort()
    points = [bb for (aa,bb) in points]
    return points
#-----------------------------------------------#
def rstlimits(drst):
    points = sorted_points(drst,hess=False)
    if len(points) == 0:
        return None, None, None, None, None, None
    # first point
    sbw, Ebw = drst[points[ 0]][0:2]
    # last point
    sfw, Efw = drst[points[-1]][0:2]
    return points[0],points[-1],sbw,sfw,Ebw,Efw
#-----------------------------------------------#
def cubic2float(cubic):
    if type(cubic) == type(""):
       cubic = cubic.lower()
    # convert to d3
    try:
        d3 = float(cubic)
    except:
        if   cubic in ["no",False,None]: d3 = None
        elif cubic in ["yes",True]     : d3 = STEP_CUBIC
        else                           : d3 = STEP_CUBIC
    # return d3
    return d3
#-----------------------------------------------#
def sd_getv0(gms=None):
    '''
    Convert mass-scaled gradient (gms) to v^(0)
    Returns None if no gradient given
    '''
    if gms is None: return None
    norm = np.linalg.norm(gms)
    v0   = [-gi/norm for gi in gms]
    return np.array(v0)
#-----------------------------------------------#
def sd_getv1(Fms,gms,v0):
    '''
    Convert mass-scaled hessian (Fms) to v^(1)
    Mass-scaled gradient and v^(0) also needed.
    Returns None if no hessian given
    '''
    if Fms is None: return None
    v0Fv0 = float( np.matrix(v0) * Fms * np.matrix(v0).transpose() )
    component_A = np.array(Fms * np.matrix(v0).transpose()).transpose()[0]
    component_B = v0Fv0*v0
    v1 = (component_A - component_B) / np.linalg.norm(gms)
    return v1
#-----------------------------------------------#
def sd_getv0v1(gms,Fms):
    '''
    Convert mass-scaled gradient and hessian to
    v^(0) and v^(1) vectors.
    Use sd_getv0 and sd_getv1 functions
    '''
    v0 = sd_getv0(gms)
    v1 = sd_getv1(Fms,gms,v0)
    return v0, v1
#-----------------------------------------------#
def sd_taylor(xms,ds,gms=None,Fms=None,v0=None,v1=None):
    '''
    Generate next structure of steepest descent path
    using Taylor series. Requires gradient in mass-scaled.
    Optionally, hessian can be given for better results.
    x(s0+ds) = x(s0) + v0 * ds + 0.5 * v1 * ds^2 + ...
    '''
    # gms,Fms --> v0, v1
    if v0 is None:
        v0 = sd_getv0(gms)
    if v1 is None:
        v1 = sd_getv1(Fms,gms,v0)
    # apply v0 and v1
    xms_new = [xi+ds*vi for xi,vi in zip(xms,v0)]
    if v1 is not None:
       xms_new = [xi+0.5*(ds**2)*vi for xi,vi in zip(xms_new,v1)]
    # return new geom
    return xms_new
#===============================================#


#===============================================#
# Functions for Page-McIver (pm) path           #
#===============================================#
def pm_dsoverdt(t, gp, lamb):
    '''
    returns ds/dt
    '''
    nelements = lamb.shape[0]
    dsdt = 0.0
    for i in range(nelements):
        dsdt += (  gp[i,0] * fncs.exp128( -lamb[i,i] * t)  )**2
    dsdt = np.sqrt( dsdt)
    return dsdt
#-----------------------------------------------#
def pm_f4root(t0,ds,gp,lamb,imethod="quad"):
    '''
    Get integral (I) and calculate:
        g(t) = delta(s) - I
    '''
    if imethod == "trap":
        # trapezoidal integration in [0,t0]
        if   t0 < 10**2: nsteps = 10**3
        elif t0 < 10**3: nsteps = 10**4
        else           : nsteps = 10**5
        dt        = 1.0 * (t0-0)/nsteps
        tvalues   = [0.0+i*dt for i in range(nsteps)]
        integrand = [pm_dsoverdt(t,gp,lamb) for t in tvalues]
        integral  = sum([y*dt for y in integrand])
    if imethod == "quad":
        integral = quad(pm_dsoverdt, 0.0, t0, args=(gp, lamb), full_output=0, limit=50)[0]
    gt = ds - integral
    return gt
#-----------------------------------------------#
def pm_gett(ds,gp,alpha,tmethod="quad"):
    '''
    returns value of t for page-mciver method
    If not able to get it, None is returned.
    alpha: eigenvalues of Fms (in diagonal matrix)
    gp   : U^T gms, with U being eigenvectors of Fms
    '''
    # (Try to) get value of t
    t0 = ds
    # find interval for root [0,t0]
    try:
        while True:
           gt = pm_f4root(t0,ds,gp,alpha,tmethod)
           if gt < 0.0: break
           t0 *= 10.0
        t = brenth(pm_f4root, t0/10.0, t0, args=(ds,gp,alpha,tmethod))
    # A problem occurs! Calculation of t value is omitted
    except:
        t = None
    return t
#-----------------------------------------------#
def pm_nextxms(xms,ds,gms,Fms,t0=None):
    '''
    Calculate next geometry using page-mciver method
    If fail in calculating t parameter, quadratic Taylor
    expasion is used.
    '''
    tmethod = "quad" # ("quad" or "trap")
    natoms = fncs.howmanyatoms(xms)
    # to numpy arrays
    if xms is not None: xms = np.matrix(np.array([xms]).transpose())
    if gms is not None: gms = np.matrix(np.array([gms]).transpose())
    if Fms is not None: Fms = np.matrix(np.array(Fms))
    # Get Hessian eigenvalues and eigenvectors
    alpha, U = np.linalg.eigh(Fms)
    # Get projection of gradient
    gp = U.transpose() * gms
    # Put eigenvectors in diagonal matrix
    alpha = np.diag(alpha)
    # check dimensions
    if gms.shape != (3*natoms,       1): exit("Problems in 'pm_nextxms'")
    if gp.shape  != (3*natoms,       1): exit("Problems in 'pm_nextxms'")
    if U.shape   != (3*natoms,3*natoms): exit("Problems in 'pm_nextxms'")
    # get t
    if t0 is None: t = pm_gett(ds,gp,alpha,tmethod)
    else         : t = t0
    # apply taylor?
    if t is None:
       xms = np.array(xms.transpose().tolist()[0])
       gms = np.array(gms.transpose().tolist()[0])
       Fms = np.array(Fms.tolist()               )
       return sd_taylor(xms,ds,gms,Fms), None
    # Get Mn matrix
    M = np.identity(3*natoms)
    for i in range(3*natoms):
        M[i,i] = (fncs.exp128(- alpha[i,i] * t) - 1.0) / alpha[i,i]
    # Get D(t)
    D = U * M * U.transpose()
    # Get dx vector
    dx_vector = D * gms
    # to normal array
    xms       = np.array(      xms.transpose().tolist()[0])
    dx_vector = np.array(dx_vector.transpose().tolist()[0])
    # Get new geom
    xms_next = xms + dx_vector
    return xms_next, t
#===============================================#





#===============================================#
# Functions related to first step of MEP        #
#===============================================#
def get_imagfreq(freqs,evecs):
    tfreq, tLms = None, None
    smallest = 0.0
    for freq,Lms in zip(freqs,evecs):
        if freq >= 0.0: continue
        if freq < smallest:
            smallest = freq
            tfreqs   = freq
            tLms     = Lms
    return tfreq, tLms
#-----------------------------------------------#
def correctdir_v0(xms,v0,idir,masses,mu):
    # Check direction of v0
    ic, sign = idir
    changedir = not intl.ics_correctdir(xms,v0,ic,sign,masses,mu)
    if changedir: v0 = -v0
    return v0
#-----------------------------------------------#
def calculate_v1(xms,v0,Fms,symbols,masses,mu,d3,spc_fnc,spc_args,parallel=False):
    natoms = fncs.howmanyatoms(xms)
    # Calculation of hessian at close structures
    xbw_3rd = fncs.ms2cc_x(sd_taylor(xms,d3,v0=-v0),masses,mu)
    xfw_3rd = fncs.ms2cc_x(sd_taylor(xms,d3,v0=+v0),masses,mu)
    t3rd_bw = (xbw_3rd,symbols,True,"thirdBw",spc_args)
    t3rd_fw = (xfw_3rd,symbols,True,"thirdFw",spc_args)
    # calculation with software
    if fncs.do_parallel(parallel):
       import multiprocessing
       pool  = multiprocessing.Pool()
       out   = [pool.apply_async(spc_fnc,args=args) for args in [t3rd_bw,t3rd_fw]]
       #outbw, outfw = Parallel(n_jobs=-1)(delayed(spc_fnc)(*inp) for inp in [t3rd_bw,t3rd_fw])
       outbw = out[0].get()
       outfw = out[1].get()
       # clean up pool
       pool.close()
       pool.join()
    else:
       outbw = spc_fnc(*t3rd_bw)
       outfw = spc_fnc(*t3rd_fw)
    xcc_bw, atnums, ch, mtp, V0_bw, gcc_bw, Fcc_bw, dummy = outbw
    xcc_fw, atnums, ch, mtp, V0_fw, gcc_fw, Fcc_fw, dummy = outfw
    # In mass-scaled
    Fms_bw = fncs.cc2ms_F(fncs.lowt2matrix(Fcc_bw),masses,mu)
    Fms_fw = fncs.cc2ms_F(fncs.lowt2matrix(Fcc_fw),masses,mu)
    # Convert numpy matrices
    Fms    = np.matrix(Fms)
    Fms_bw = np.matrix(Fms_bw)
    Fms_fw = np.matrix(Fms_fw)
    # Matriz of third derivatives
    m3der = (Fms_fw - Fms_bw) / (2.0*d3)
    # Calculation of curvature, i.e. v^(1) (or small c)  A = B * c --> c = B^-1 * A
    LF = np.matrix([v0]).transpose()
    A  = m3der*LF - float(LF.transpose()*m3der*LF) * LF
    B  = 2.0 * float( LF.transpose() * Fms * LF) * np.identity(3*natoms) - Fms
    v1 = np.linalg.inv(B) * A
    # Convert to (normal) array
    v1 = np.array(v1.transpose().tolist()[0])
    return v1
#-----------------------------------------------#
def mep_first(xcc,gcc,Fcc,symbols,masses,var_first,spc_fnc,spc_args,drst={},parallel=False):

    #global PARALLEL, delayed, multiprocessing, Parallel
    #PARALLEL, delayed, multiprocessing, Parallel = fncs.set_parallel(parallel)
    ds,mu,cubic,idir = var_first

    # convert cubic variable to d3 (in bohr)
    d3 = cubic2float(cubic)

    # mass-scaled
    xcc = fncs.shift2com(xcc,masses)
    xms = fncs.cc2ms_x(xcc,masses,mu)
    gms = fncs.cc2ms_g(gcc,masses,mu)
    Fms = fncs.cc2ms_F(Fcc,masses,mu)

    # Data in backup?
    if TSLABEL in drst.keys():
        si, E, xms2, gms, Fms, v0, v1, t = drst[TSLABEL]
        gms = np.array(gms)
        Fms = np.array(Fms)
        v0  = np.array(v0 )
        if v1 is not None: v1  = np.array(v1 )
        same = fncs.same_geom(xms,xms2)
        if not same: raise Exc.RstDiffTS(Exception)
    else:
        # Calculation of v0
        freqs, evals, evecs = fncs.calc_ccfreqs(Fcc,masses,xcc)
        ifreq, Lms          = get_imagfreq(freqs,evecs)
        v0  = fncs.normalize_vec(Lms)
        v0  = correctdir_v0(xms,v0,idir,masses,mu)

        # Calculation of v1
        if d3 is None: v1 = None
        else         : v1 = calculate_v1(xms,v0,Fms,symbols,masses,mu,d3,spc_fnc,spc_args,parallel)
    # Final structures
    if d3 is None:
       xms_bw = sd_taylor(xms,ds,v0=-v0)
       xms_fw = sd_taylor(xms,ds,v0=+v0)
    else:
       xms_bw = sd_taylor(xms,ds,v0=-v0,v1=-v1)
       xms_fw = sd_taylor(xms,ds,v0=+v0,v1=+v1)

    return (xms, gms, Fms), (v0, v1), (xms_bw, xms_fw)
#===============================================#


#===============================================#
# Functions related to the calculation of the   #
# path from a point != from saddle point        #
#===============================================#
def steepest(xcc0,s0,symbols,masses,pathdata,spc_fnc,spc_args,drst={},lFms=None,pf="%i"):
    '''
    x0 should not be a saddle point
    pathdata: tuple made of:
        * method
          pm for Page-McIver
          es for Euler
        * mu
        * ds
        * nsteps
        * hsteps
        * epse
        * epsg
    pf: point format
    '''
    # expand pathdata
    method,mu,ds,sfinal,hsteps,epse,epsg = pathdata
    nsteps = int(round(abs(sfinal/ds)))
    # Check some staff
    if lFms is None and method == "pm": raise Exc.SDpmNoHess(Exception)
    # loop
    s_i   = s0
    xcc_i = xcc0
    listE = []
    listg = []
    for step in range(1,nsteps+1):
        point = pf%step
        # calculate hessian?
        if step%hsteps == 0: bhessian = True
        else               : bhessian = False
        # previously calculated vs single-point calculation
        if point in drst.keys():
            s_ii, E_i, xms_i, gms_i, Fms_i, v0, v1, t_i = drst[point]
            # in cc
            gcc_i = fncs.ms2cc_g(gms_i,masses,mu)
            # compare geometry
            exception = Exc.RstDiffPoint(Exception)
            exception._var = (point,s_i)
            if abs(s_ii-s_i) > EPS_MEPS: raise exception
            same = fncs.same_geom(fncs.cc2ms_x(xcc_i,masses,mu),xms_i)
            if not same: raise exception
        else:
            spc_data = spc_fnc(xcc_i,symbols,bhessian,point,spc_args)
            xcc_i, atnums, ch, mtp, E_i, gcc_i, Fcc_i, dummy = spc_data
            # Data in mass-scaled and correct format
            Fcc_i = fncs.lowt2matrix(Fcc_i)
            xms_i = fncs.cc2ms_x(xcc_i,masses,mu)
            gms_i = fncs.cc2ms_g(gcc_i,masses,mu)
            Fms_i = fncs.cc2ms_F(Fcc_i,masses,mu)
            t_i   = None
        # keep data for eps
        listE.append(E_i)
        listg.append(fncs.norm(gcc_i))
        # Last Fms
        if Fms_i not in [[],None]: lFms = Fms_i
        # Calculate new geometry (step+1)
        if   method == "es":
           xms_j, t_i = sd_taylor(xms_i,ds,gms=gms_i), None # to check format
        elif method == "pm":
           xms_j, t_i = pm_nextxms(xms_i,ds,gms=gms_i,Fms=lFms,t0=t_i)
        # check real ds value
        diff_ds2 = sum([(ii-jj)**2 for (ii,jj) in zip(xms_j,xms_i)]) - ds**2
        if diff_ds2 > EPS_FLOAT: raise Exc.SDdiffds(Exception)
        # Yield updated data for step
        yield (point,s_i), E_i, xms_i, gms_i, Fms_i, t_i
        # update
        xcc_i = fncs.ms2cc_x(xms_j,masses,mu)
        if   s_i >= 0.0: s_i += ds
        else           : s_i -= ds
        # Check eps
        if len(listE) > 2 and step > 2*hsteps:
            diffE = abs(listE[-1]-listE[-2])
            #diffg = abs(listg[-1]-listg[-2])
            crit1 = diffE     < epse
            crit2 = listg[-1] < epsg
            if crit1 and crit2: break
#===============================================#
