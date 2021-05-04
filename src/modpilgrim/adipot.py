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
| Sub-module :  adipot             |
| Last Update:  2021/04/20 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#=======================================================#
import common.fncs           as fncs
import common.Exceptions     as Exc
#-------------------------------------------------------#
import modpilgrim.steepdesc  as sd
#-------------------------------------------------------#
from   common.criteria   import EPS_CCICZPE
from   common.criteria   import EPS_CCIC
from   common.Molecule   import Molecule
from   common.physcons   import KCALMOL
from   common.Spline     import VadiSpline
#=======================================================#

#=======================================================#
def ccvsic_checknfreqs(lcc_frqs,lic_frqs):
    if len(lic_frqs) == 0: return True
    for ccfrqs,icfrqs in zip(lcc_frqs,lic_frqs):
        if len(ccfrqs) != len(icfrqs): return False
    return True
#-------------------------------------------------------#
def ccvsic_checkts(data_x,lcc_frqs,lic_frqs):
    # if no ic-freqs, nothing to compare
    if len(lic_frqs) == 0: return True
    # get idx of TS
    idxTS = None
    for idx,s_i in enumerate(data_x):
        if s_i == 0.0: idxTS=idx; break
    # TS not found?
    if idxTS is None: raise Exception
    # Compare date
    ccfrqs = lcc_frqs[idxTS]
    icfrqs = lic_frqs[idxTS]
    # different length?
    if len(ccfrqs) != len(icfrqs): return False
    # maximum frequency difference (in 1/cm)
    freq_diff = max([abs(f1-f2) for f1,f2 in zip(ccfrqs,icfrqs)])
    freq_diff = fncs.afreq2cm(freq_diff)
    if freq_diff > EPS_CCIC: return False
    return True
#-------------------------------------------------------#
def ccvsic_checkzpe(lcc_tzpe,lic_tzpe):
    if len(lic_tzpe) == 0: return True
    for zpe_cc,zpe_ic in zip(lcc_tzpe,lic_tzpe):
        diff = (zpe_ic-zpe_cc)*KCALMOL
        if diff < -EPS_CCICZPE: return False
    return True
#-------------------------------------------------------#
def path2vadi(tcommon,drst,Eref=None,ics=None,boolint=False,lowfq={},freqscal=1.0,icsbw=None,icsfw=None,symmetry=None):
    '''
    lowfq: in case of imaginary frequencies
    '''
    # boolint as boolean
    if   boolint in [False,"no","No","n","N",None]: boolint = False
    elif boolint in [True,"yes","YES","y","Y"]    : boolint = True
    # check there are internal coordinates if required
    if boolint and ics in [None,False,[]]: raise Exc.NoICS(Exception)
    # expand data in tcommon
    ch,mtp,atnums,masses,mu = tcommon
    # Sorted labels (by s)
    slabels = sd.sorted_points(drst,hess=True)
    # Reference energy
    if Eref is None: lbw, lfw, sbw, sfw, Eref, Efw = sd.rstlimits(drst)
    # Independent variable
    data_x = [drst[label][0] for label in slabels]
    # mep energy (relative value)
    listV0 = [drst[label][1]-Eref for label in slabels]

    # Initializing data
    lcc_tzpe, lcc_frqs, lcc_Vadi = [], [], []
    lic_tzpe, lic_frqs, lic_Vadi = [], [], []
    dMols = {}

    # Updating data
    for label in slabels:
        # data in drst
        s_i, E_i, xms_i,gms_i,Fms_i,v0_i,v1_i,t_i = drst[label]
        # project gradient
        if s_i == 0.0: bool_pg = False
        else         : bool_pg = True
        # lowfq
        if   s_i == 0.0: dlowfq = {}
        elif s_i  < 0.0: dlowfq = lowfq.get("bw",{})
        elif s_i  > 0.0: dlowfq = lowfq.get("fw",{})
        # mass-scaled --> Cartesian coords
        xcc = fncs.ms2cc_x(xms_i,masses,mu)
        gcc = fncs.ms2cc_g(gms_i,masses,mu)
        Fcc = fncs.ms2cc_F(Fms_i,masses,mu)
        # create Molecule instance
        mol = Molecule()
        mol.setvar(xcc=xcc,gcc=gcc,Fcc=Fcc)
        mol.setvar(atonums=atnums,masses=masses)
        mol.setvar(ch=ch,mtp=mtp,V0=E_i)
        mol.setvar(fscal=freqscal)
        if symmetry is not None:
            mol.setvar(pgroup=symmetry[0])
            mol.setvar(rotsigma=symmetry[1])
        mol.prepare()
        mol.setup(mu,projgrad=bool_pg)
        mol.clean_freqs("cc")       # it may be needed with orca
        mol.deal_lowfq(dlowfq,"cc") # deal with low frequencies
        mol.ana_freqs("cc")         # calculate zpe
        # append data
        lcc_tzpe.append(float(mol._cczpe)   )
        lcc_Vadi.append(mol._ccV1 - Eref    )
        lcc_frqs.append(list(mol._ccfreqs))
        # internal coordinates
        if boolint:
           if   s_i < 0.0 and icsbw is not None: mol.icfreqs(icsbw,bool_pg)
           elif s_i > 0.0 and icsfw is not None: mol.icfreqs(icsfw,bool_pg)
           else                                : mol.icfreqs(ics  ,bool_pg)
           mol.deal_lowfq(dlowfq,"ic")
           mol.ana_freqs("ic")
           # append data
           lic_tzpe.append(float(mol._iczpe) )
           lic_Vadi.append(mol._icV1 - Eref  )
           lic_frqs.append(list(mol._icfreqs))
        # save instance
        dMols[label] = (s_i,mol)

    # package data
    tuple_cc = (data_x,lcc_frqs,lcc_tzpe)
    tuple_ic = (data_x,lic_frqs,lic_tzpe)
    # Generate splines
    Vadi_cc = VadiSpline(data_x,lcc_Vadi)
    if boolint: Vadi_ic = VadiSpline(data_x,lic_Vadi)
    else      : Vadi_ic = None
    # Return data
    return dMols, Vadi_cc, Vadi_ic, tuple_cc, tuple_ic, listV0
#=======================================================#


