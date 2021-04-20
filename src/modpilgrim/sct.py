'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.2
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
| Sub-module :  sct                |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains functions for the
calculation of the SCT correction factor

By default, v^(1) vectors are calculated through
the hessian matrix. However, they can be calculated
using numerical differences in the gradients.
To do so, use the function 'get_numv1' to get
the dictionary 'dv1'. This dictionary should be
used in 'get_sct_part1' and 'get_sct' functions.
'''

#===============================================#
import numpy        as np
from   scipy.integrate import fixed_quad
from   scipy.integrate import simps
#-----------------------------------------------#
import modpilgrim.steepdesc    as sd
#-----------------------------------------------#
import common.fncs         as fncs
import common.interpolate  as intrpl
import common.physcons     as pc
#-----------------------------------------------#
from   common.Spline          import Spline
from   common.criteria        import EPS_MEPE
#===============================================#

INTRPLMODE = "linear" # linear/cubic
INTRPLNUM  = 0        # 0/1/2...

#===============================================#
#     Numeric calculation of v^(1) vectors      #
#===============================================#
def get_numv1(drst):
    # list of point
    pts1 = sd.sorted_points(drst,hess=False)
    pts2 = sd.sorted_points(drst,hess=True )
    # get some lists from drst
    pts1  = sd.sorted_points(drst,hess=False)
    pts2  = sd.sorted_points(drst,hess=True )
    svals = [drst[label][0] for label in pts1]
    grads = [drst[label][3] for label in pts1]
    # positions in pts1 of points in pts2
    hess_idxs   = [idx for idx,label in enumerate(pts1) if label in pts2 ]
    hess_idxsBW = [idx for idx,label in enumerate(pts1) if label in pts2 and drst[label][0]<0.0]
    hess_idxsFW = [idx for idx,label in enumerate(pts1) if label in pts2 and drst[label][0]>0.0]
    # calculate list of numeric v^(1) vectors
    if False: dv1 = numv1_calcA(pts1,svals,grads,hess_idxsBW,hess_idxsFW)
    if True : dv1 = numv1_calcB(pts1,svals,grads,hess_idxs) # faster than numv1_calcA
    return dv1
#-----------------------------------------------#
def numv1_calcA(all_points,all_svals,grads,hess_idxsBW,hess_idxsFW):
    '''
    svals: the list of s values
    grads: the list of gradients
    hess_idxs: a list of indices indicating where v^(1)
               has to be calculated
               this is given in two list, one for the bw and
               another for the fw part of the MEP
               why? although each coordinate of the gradient is continuous along the mep,
               this behavior is not found for the v^(0) vector.
    '''
    dv1 = {}
    for indices in [hess_idxsBW,hess_idxsFW]:
        for idx in indices:
          # label for this point
          label = all_points[idx]

          # Choose the points
          if   idx == indices[ 0]: idx_a, idx_b, idx_c = idx  , idx+1, idx+2
          elif idx == indices[-1]: idx_a, idx_b, idx_c = idx-2, idx-1, idx  
          else                   : idx_a, idx_b, idx_c = idx-1, idx  , idx+1

          # s for the three points of mep
          si = all_svals[idx]
          sa = all_svals[idx_a]
          sb = all_svals[idx_b]
          sc = all_svals[idx_c]

          # Get v^(0) at each point
          v0_a = - fncs.normalize_vec(grads[idx_a])
          v0_b = - fncs.normalize_vec(grads[idx_b])
          v0_c = - fncs.normalize_vec(grads[idx_c])

          # Parabolic function for each coordinate
          ncoords = len(v0_a)
          v1      = []
          for coord in range(ncoords):
              x1, x2, x3 = sa, sb, sc
              y1, y2, y3 = v0_a[coord], v0_b[coord], v0_c[coord]
              # Get a in a*x^2 + b*x + c
              num_a = (y1-y2) * (x1-x3) - (y1-y3) * (x1-x2)
              den_a = (x1**2 - x2**2) * (x1-x3) - (x1**2 - x3**2) * (x1-x2)
              a = num_a / den_a
              # Get b in b*x^2 + b*x + c
              b = ((y1 - y2) - a * (x1**2 - x2**2)) / (x1-x2)
              # Get c in b*x^2 + b*x + c
              c = y1 - a*x1**2 - b *x1
              # Get derivative at x=s (y' = 2*a*x+b)
              v1_coord = 2 * a * si + b
              v1.append(v1_coord)
          dv1[label] = np.array(v1)
    return dv1
#-----------------------------------------------#
def numv1_calcB(points,svals,grads,hess_idxs):
    '''
    average of linear interpolation
    it is equivalent to numv1_calcA, but faster
    _cu for current point
    _m1 for previous point to current
    _m2 for previous point to _m1
    _p1 for next point with regard to current
    _p2 for next point with regard to _p1
    et cetera
    '''
    dv1 = {}
    saddle_idx = None

    for idx,label in enumerate(points):
        if idx not in hess_idxs: continue
        # Current point
        s_cu  = svals[idx]
        g_cu  = np.array(grads[idx])
        # Skip if saddle point
        if s_cu == 0.0: continue
        # v0 vector
        v0_cu = - fncs.normalize_vec(g_cu)
        # Previous point
        if idx != 0:
          s_m1  = svals[idx-1]
          g_m1  = np.array(grads[idx-1])
          v0_m1 = - fncs.normalize_vec(g_m1)
          # Get left derivative
          v1_left = (v0_cu - v0_m1) / (s_cu-s_m1)
        else:
          v1_left = None
        # Next point
        if idx != len(points) -1:
          s_p1  = svals[idx+1]
          g_p1  = np.array(grads[idx+1])
          v0_p1 = - fncs.normalize_vec(g_p1)
          # Get right derivative
          v1_right = (v0_p1 - v0_cu) / (s_p1-s_cu)
        else:
          v1_right = None
        # average
        if   v1_left  is None: v1 = v1_right
        elif v1_right is None: v1 = v1_left
        else:                  v1 = (v1_right+v1_left)/2.0
        # save
        dv1[label] = v1
    return dv1
   ## Interpolate v1 for saddle point
   #if saddle_idx is not None:
   #   idx2 = hess_idxs.index(saddle_idx)
   #   exit()

   #   s_m2 = svals[saddle_idx-2]
   #   s_m1 = svals[saddle_idx-1]
   #   s_cu = svals[saddle_idx+0]
   #   s_p1 = svals[saddle_idx+1]
   #   s_p2 = svals[saddle_idx+2]

   #   v1_m2 = dv1[points[saddle_idx-2]]
   #   v1_p2 = dv1[points[saddle_idx+2]]
   #   exit()

   #   v1_m1 = (s_m1-s_m2)/(s_p2-s_m2) * (v1_p2 - v1_m2) + v1_m2
   #   v1_cu = (s_cu-s_m2)/(s_p2-s_m2) * (v1_p2 - v1_m2) + v1_m2
   #   v1_p1 = (s_p1-s_m2)/(s_p2-s_m2) * (v1_p2 - v1_m2) + v1_m2

   #   dv1[points[saddle_idx-1]] = v1_m1
   #   dv1[points[saddle_idx+0]] = v1_cu
   #   dv1[points[saddle_idx+1]] = v1_p1
#===============================================#

def get_dtbar_ds(s_list, tbar_list):
    '''
    OLD FUNCTION TO CALCULATE dt/ds
    NO LONGER USED!
    --------------------------------
    Function to obtain the numerical derivative
    of tbar with regard to the MEP coordinate s
    PS: tbar is the classical turning point associated
        to the whole set of harmonic modes
    '''
    num_points = len(s_list)
    dtbards_list = []

    for idx in range(num_points):
        s     = s_list[idx]
        # Get dt/ds
        if idx == 0:
            s1, s2, s3 = s_list[0:3]
            t1, t2, t3 = tbar_list[0:3]
            dtbar_ds = (-3*t1 + 4*t2 - t3) / (s3-s1)
        elif idx == num_points-1:
            sm2, sm1, sm = s_list[idx-2:idx+1]
            tm2, tm1, tm = tbar_list[idx-2:idx+1]
            dtbar_ds = (tm2 - 4*tm1 + 3*tm) / (sm-sm2)
        else:
            s_previous = s_list[idx-1]
            s_next     = s_list[idx+1]
            t_previous = tbar_list[idx-1]
            t_next     = tbar_list[idx+1]
            dtbar_ds = (t_next - t_previous)/(s_next-s_previous)
        # Append data
        dtbards_list.append(dtbar_ds)
    return dtbards_list


#===============================================#
#   Calculation of theta and kappa integrand    #
#===============================================#
def theta_integrand(s_i,E,svals,lmueff,VadiSpl):
    '''
    lmueff: a float (for ZCT) or a list (for SCT)
    '''
    #-------------------------------------------------#
    def theta_integrand_float(s_i,E,svals,lmueff,VadiSpl):
        if type(lmueff) == float:
           mu_eff = lmueff
        else:
           mu_eff = max(0.0,intrpl.interpolate(svals,lmueff,s_i))
        return np.sqrt( 2*mu_eff* abs(E-VadiSpl(s_i)) )
    #-------------------------------------------------#
    # treat as list
    if type(s_i) in (float,np.float64,np.float128):
        return  theta_integrand_float(s_i,E,svals,lmueff,VadiSpl)
    else:
        return [theta_integrand_float(x_i,E,svals,lmueff,VadiSpl) for x_i in s_i]
#-----------------------------------------------#
def get_theta(E,svals,lmueff,VadiSpl):
    # at the top?
    sAG, vAG = VadiSpl.get_max()
    if abs(E-vAG) < EPS_MEPE: return 0.0, [], 0.0
    # return points
    rpoints = VadiSpl.returnpoints(E)
    # calculate theta
    theta1, theta2 = 0.0, 0.0
    for (si,sj) in rpoints:
        args = (E,svals,lmueff,VadiSpl)
        # integrate with two methods
        integral1 = fncs.intg_gau( theta_integrand, si, sj, n= 80, args=args)
        integral2 = fncs.intg_trap(theta_integrand, si, sj, n=160, args=args)
        # add to theta
        theta1 += integral1
        theta2 += integral2
    # choose theta from integrations
    theta = min(theta1,theta2)
    if theta != 0.0: diff = 100*abs(theta1-theta2)/theta
    else           : diff = 0.0
  # # if difference is significant, I trust the trapezoidal integration the most
  # if diff > 10.0: theta = theta2
    # return
    return theta, rpoints, diff
#-----------------------------------------------#
def kappa_int(VAG,E,theta,beta):
    pE = 1.0 / (1.0 + fncs.exp128(2.0*theta))
    kappa_integrand = pE * np.sinh(beta*(VAG-E))
    return kappa_integrand
#-----------------------------------------------#
def gauquad_pointsweights(ntp,x0=-1.0,xn=1.0):
    points, weights = np.polynomial.legendre.leggauss(ntp)
    suma    = (xn+x0)/2.0
    resta   = (xn-x0)/2.0
    points  = [resta*xi+suma for xi in points]
    weights = [wi*resta for wi in weights]
    return points, weights
#-----------------------------------------------#
def kappa_integral1(E_list,probs,weights,beta,VAG):
    '''
    integral from E0 to VAG
    '''
    integral = 0.0
    INTGRND  = []
    for idx in range(len(E_list)):
        E_i       = E_list[idx]
        w_i       = weights[idx]
        pE        = probs[idx]
        integrand =  pE * fncs.exp128(-beta*(E_i-VAG)) * beta
        INTGRND.append(integrand)
        integral += w_i*integrand   
        del integrand
    return integral, INTGRND
#-----------------------------------------------#
def kappa_integral2(E_list,probs,weights,beta,VAG):
    integral = 0.0
    INTGRNDX = []
    INTGRNDY = []
    for idx in range(len(E_list)):
        E_i       = E_list[idx]
        w_i       = weights[idx]
        pE        = probs[idx]
        integrand =  (1.0-pE) * fncs.exp128(-beta*(2*VAG-E_i-VAG)) * beta
        INTGRNDX.append(2*VAG-E_i)
        INTGRNDY.append(integrand)
        integral += w_i*integrand   
        del integrand
    return integral, INTGRNDX[::-1], INTGRNDY[::-1]
#-----------------------------------------------#
def kappa_integral3(E0,VAG,beta):
    return fncs.exp128(-beta*(2*VAG-E0-VAG))
#===============================================#




#===============================================#
def get_sct_part1(points,VadiSpl,E0=None):
    # Check E0 value
    s_bw, E0_bw = VadiSpl.get_alpha()
    s_fw, E0_fw = VadiSpl.get_omega()
    if E0 is None: E0 = max(E0_bw,E0_fw)
    return E0
#-----------------------------------------------#
def get_sct_part2(dMols,points,dv1={},case="cc",INTRPL=(INTRPLMODE,INTRPLNUM)):
    '''
    calculates BmF and tbar
    at each point of the path
    if dv1 == {}, analytic v1 vectors will be used
    calculates effective mass
    '''

    intrplmode, intrplnum = INTRPL
    if intrplmode not in ["linear","cubic"]: intrplmode = INTRPLMODE
    if type(intrplnum) != int: intrplnum = INTRPLNUM
    if      intrplnum   < 0  : intrplnum = INTRPLNUM
    if      intrplnum   > 2  : intrplnum = INTRPLNUM

    svals = [dMols[point][0] for point in points]

    #------------------#
    # Points to ignore #
    #------------------#
    # localize transition state
    idxts = None
    for idx,label in enumerate(points):
        s_i = dMols[label][0]
        if s_i == 0.0: idxts = idx
    # ignore points
    toignore = []
    if idxts is not None:
       # ignore TS
       toignore += [idxts]
       # ignore closest points to the TS
       npbw = len(points[0:idxts] )
       npfw = len(points[idxts+1:])
       if npbw - intrplnum >= 2: toignore += [idxts-ii  for ii in range(1,intrplnum+1,1)]
       if npfw - intrplnum >= 2: toignore += [idxts+ii  for ii in range(1,intrplnum+1,1)]

    #-------------------------------#
    # calculate Bmf, tbar and dtbar #
    #-------------------------------#
    lkappa, ltbar = [], []
    for idx,label in enumerate(points):
        s_i = dMols[label][0]
        mu  = dMols[label][1]._mu
        if idx in toignore:
           lkappa.append(None)
           ltbar.append(None)
           continue
        # Get v1
        if label in list(dv1.keys()): v1 = dv1[label]
        else: v1 = sd.sd_getv0v1(dMols[label][1]._gms,dMols[label][1]._Fms)[1]
        # frequencies and evectors
        if   case == "cc":
           freqs = dMols[label][1]._ccfreqs
           evecs = dMols[label][1]._ccFevecs
        elif case == "ic":
           freqs = dMols[label][1]._icfreqs
           evecs = dMols[label][1]._icFevecs
        # Calculate Bmf's
        bmfs    = [ - fncs.sign(s_i) * np.dot(Lm,v1) for Lm in evecs]
        # Calculate turning points
        turnpts = [fncs.afreq2turnpoint(freq,mu) for freq in freqs]
        # Calculate kappa
        kappa = np.sqrt(sum([bmf**2 for bmf in bmfs]))
        # Calculate tbar
        tbar = sum( [ (bmf/tp/tp)**2 for (bmf,tp) in zip(bmfs,turnpts)] ) ** (-0.25) * np.sqrt(kappa)
        # Append data
        lkappa.append( kappa )
        ltbar.append( tbar )

    # Interpolate Nones
    lkappa = intrpl.interpolate_nones(svals,lkappa,mode=intrplmode)
    ltbar  = intrpl.interpolate_nones(svals,ltbar ,mode=intrplmode)

    # Derivative via spline
    ldtbar = [intrpl.interpolate(svals,ltbar,s_i,d=1)[1] for s_i in svals]

    #-------------------------------#
    # calculation of effective mass #
    #-------------------------------#
    lfs   = [] # f(s)
    tsidx = None
    for idx,s_i in enumerate(svals): 
        kappa    = lkappa[idx]
        tbar     = ltbar[idx]
        dtbar    = ldtbar[idx]
        if s_i == 0.0: tsidx = idx
        # Calculare f(s) according to eq (14) - JAmChemSoc(1993)_115_2408
        if idx in toignore:
           lfs.append(None)
        else:
           exparg = -2*kappa*tbar - (kappa*tbar)**2 + dtbar**2
           fs     = min(fncs.exp128(exparg),1.0)
           lfs.append(fs)
    # Interpolate Nones
    lfs = intrpl.interpolate_nones(svals,lfs,mode=intrplmode)
    # set to zero possible negative values due to interpolation
    for idx in toignore: lfs[idx] = max(0.0,lfs[idx])
    # effective mu in au
    lmueff = [mu*fs for fs in lfs]

    return svals, lkappa, ltbar, ldtbar, mu, lmueff, toignore
#-----------------------------------------------#
def get_sct_part3(svals,lmueff,VadiSpl,afreq,lEquant,E0,VAG,temps):
    '''
    discrete sum from E0 to VAG
    Surface Science, 164, 558-588 (1985)
    '''
    lCOEF = []
    # Only energies between E0 and VAG
    lEquant = [E for E in lEquant if  E0-EPS_MEPE <= E <= VAG+EPS_MEPE]
    # calculate theta and p(E) for each energies
    thetuple    = [get_theta(En,svals,lmueff,VadiSpl) for En in lEquant]
    thetas      = [theta for theta,rps,diff in thetuple]
    rpoints     = [rps   for theta,rps,diff in thetuple]
    diffs       = [diff  for theta,rps,diff in thetuple]
    # Probability
    probs  = [1.0 / (1.0 + fncs.exp128(2.0*theta)) for theta in thetas]
    # calculate coefficient
    dEndn = pc.HBAR * afreq
    for idx,T in enumerate(temps):
       beta = 1.0 / (pc.KB*T)
       kappa = 0.0
       for En,pE in zip(lEquant,probs):
           kappa += dEndn * pE * fncs.exp128(-beta*(En-VAG)) * beta
       lCOEF.append(kappa)
    return lCOEF, lEquant, probs, rpoints, diffs
#-----------------------------------------------#
def get_sct_part4(svals,lmueff,VadiSpl,E0):
    '''
    Calculation of tunnelling probabilities
      if mueff is float (=mu) --> ZCT
      if mueff is a list      --> SCT
    '''
    sAG, VAG = VadiSpl.get_max()
    # List of energies for Gaussian quadrature
    E_list,weights = gauquad_pointsweights(80,E0,VAG)
    # Calculate theta_ZCT or theta_SCT (T-independent)
    thetuple    = [get_theta(E,svals,lmueff,VadiSpl) for E in E_list]
    thetas      = [theta for theta,rps,diff in thetuple]
    rpoints     = [rps   for theta,rps,diff in thetuple]
    diffs       = [diff  for theta,rps,diff in thetuple]
    # Probabilities
    probs       = [1.0 / (1.0 + fncs.exp128(2.0*theta)) for theta in thetas]
    # Calculate theta also at E=E0
    theta0,rps0,diff0 = get_theta(E0 ,svals,lmueff,VadiSpl)
    prob0   = 1.0 / (1.0 + fncs.exp128(2.0*theta0 ))
    # return data
    return weights, E_list, probs, rpoints, diffs, (prob0,rps0)
#-----------------------------------------------#
def get_sct_part5(E_list,probs,weights,E0,VAG,temps,discrete=None,qrc_Elim=None):
    '''
    calculates the correction factor
    discrete = a list with the values for I1
    '''
    if qrc_Elim is None: qrc_Elim = float("inf")
    lCOEFs  = []
    lRTEs   = []
    lINTGR  = []
    lIi     = []
    bqrc = [False for T in temps]
    for idx,T in enumerate(temps):
       beta = 1.0 / (pc.KB*T)
       # case (a) E0 < VAG
       if E0 < VAG:
          # Integral I1 (from E0 to VAG)
          I1, INTGRND1Y = kappa_integral1(E_list,probs,weights,beta,VAG)
          INTGRND1X = E_list
          # Integral I2 (from VAG to 2VAG-E0)
          I2, INTGRND2X, INTGRND2Y = kappa_integral2(E_list,probs,weights,beta,VAG)
          # Integral I3 (from 2VAG-E0 to infty)
          I3 = kappa_integral3(E0,VAG,beta)
          # Get representative Tunnelling Energy
          INTGRNDX = INTGRND1X+INTGRND2X
          INTGRNDY = INTGRND1Y+INTGRND2Y
          spl = Spline(INTGRNDX,INTGRNDY)
          spl.find_xtr("max")
          RTE, dummy = spl.get_max()
          # If qrc --> overwrite value
          if discrete is not None and RTE < qrc_Elim:
             I1 = discrete[idx]
             bqrc[idx] = True
          # Tunnelling coefficient
          COEF = I1+I2+I3
       # case (b) E0 >= VAG 
       else:
          COEF     = 1.0
          I1,I2,I3 = None,None,None
          RTE      = None
          INTGRNDX = None
          INTGRNDY = None
       # append data
       lCOEFs.append(COEF)
       lIi.append( (I1,I2,I3) )
       lRTEs.append(RTE)
       lINTGR.append((INTGRNDX,INTGRNDY))
    return lCOEFs, lIi, lRTEs, lINTGR, bqrc
#===============================================#
def get_sct(dMols,points,VadiSpl,temps,dv1={}):
    pass
#===============================================#
